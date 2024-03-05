library(ggplot2)
library("ggsci")
library(doSNOW)
source("mbDecoda.R")
source("ZINB_sim.R")
source("omnibus.R")

ncores=80
n.sim=100
K=100
n1=25
n2=25
sig.prob=0.2
zi.prob = rep(c(0,0.1,0.3,0.5),each=2)
bias = rep(c("small","large"),4)
sim.seed=matrix(1:(n.sim*8),n.sim,8)
confounder = F

DATA=list()
for (i in 1:8) {
  zi.prob.i=zi.prob[i]
  bias.i=bias[i]
  for (j in sim.seed[,i]) {
    set.seed(j)
    DATA[[j]]=ZINB_sim(seed=j,K=K,n1=n1,n2=n2,p=sig.prob,bias =bias.i,zi=zi.prob.i,confounder =confounder)
  }
}



cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
edgeR=foreach(i = 1:length(DATA), .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(DESeq2)
  library(edgeR)
  library(Wrench)
  
  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["confounder"]]
  
  
  #wrench+edgeR
  W <- wrench( t(count), condition=group)
  compositionalFactors <- W$ccf
  edgerobj <- DGEList( counts=t(count),
                       group = as.matrix(group),
                       norm.factors=compositionalFactors )
  design.mat=model.matrix(~ 0 + edgerobj$samples$group)
  colnames(design.mat)=levels(edgerobj$samples$group)
  
  edgerobj=estimateDisp(edgerobj, design.mat)
  suppressWarnings(out_edgeR <- try(glmQLFTest(glmQLFit(edgerobj, design.mat), contrast=c(1, -1)), 
                                    silent = TRUE))
  if (inherits(out_edgeR, "try-error")) {
    fdr_edgeR=0; power_edgeR=0
  }else{
    id_edgeR <- which(out_edgeR[["table"]][["PValue"]]<0.05)
    power_edgeR=length(intersect(id_edgeR,id))/length(id)
    fdr_edgeR=1-length(intersect(id_edgeR,id))/length(id_edgeR)
  }
  
  c(power_edgeR, fdr_edgeR)
}

stopCluster(cl)


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
omnibus=foreach(i = 1:length(DATA), .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(pscl)
  library(matrixStats)
  
  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["confounder"]]
  
  
  #omnibus
  suppressWarnings(out_omnibus <- try(ZISeq(t(count), as.data.frame(group), 
                                            size.factor = NULL,                                  # Normalization (GMPR)
                                            winsor = F, winsor.qt = 0.97,                     # Winsorization
                                            grp.name = 'group', adj.name = NULL, 
                                            method = 'omnibus',
                                            filter = F, prev.filter = 0, ct.cutoff = 10),     # Filter
                                      silent = TRUE))
  if (inherits(out_omnibus, "try-error")) {
    fdr_omnibus=0; power_omnibus=0
  }else{
    which(p.adjust(out_omnibus[["result"]][["p.value"]],"BH")<0.05)
    id_omnibus <- which(p.adjust(out_omnibus[["result"]][["p.value"]],"BH")<0.05)
    power_omnibus=length(intersect(id_omnibus,id))/length(id)
    fdr_omnibus=1-length(intersect(id_omnibus,id))/length(id_omnibus)
  }
  
  
  c(power_omnibus, fdr_omnibus)
}

stopCluster(cl)


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
DESeq2=foreach(i = 1:length(DATA), .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(DESeq2)
  library(Wrench)
  
  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["confounder"]]
  
  #wrench+DESeq2
  W <- wrench( t(count), condition=group)
  normalizationFactors <- W$nf
  deseq.obj <- try(DESeqDataSetFromMatrix(countData = t(count),
                                          DataFrame(group=factor(group)),
                                          ~ group ),silent = TRUE)
  try(DESeq2::sizeFactors(deseq.obj) <- normalizationFactors,silent = TRUE)
  suppressWarnings(out_DESeq2 <- try(results(DESeq(deseq.obj), pAdjustMethod = "BH"), 
                                     silent = TRUE))
  if (inherits(out_DESeq2, "try-error")) {
    fdr_DESeq2=0; power_DESeq2=0
  }else{
    id_DESeq2 <- which(out_DESeq2@listData[["padj"]]<=0.05)
    power_DESeq2=length(intersect(id_DESeq2,id))/length(id)
    fdr_DESeq2=1-length(intersect(id_DESeq2,id))/length(id_DESeq2)
  }
  
  
  c(power_DESeq2, fdr_DESeq2)
}

stopCluster(cl)


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
mbDecoda=foreach(i = 1:length(DATA), .combine = 'cbind') %dopar% {
  library(glmmTMB)
  
  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["confounder"]]
  
  #mbDecoda
  suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                       , silent = TRUE))
  if (inherits(out_mbDecoda, "try-error")) {
    fdr_mbDecoda=0; power_mbDecoda=0
  }else{
    id_mbDecoda=which(out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]<0.05)
    power_mbDecoda=length(intersect(id_mbDecoda,id))/length(id)
    fdr_mbDecoda=1-length(intersect(id_mbDecoda,id))/length(id_mbDecoda)
  }
  
  
  c(power_mbDecoda, fdr_mbDecoda)
}

stopCluster(cl)

simlist=rbind(mbDecoda[1,],omnibus[1,],edgeR[1,],DESeq2[1,],mbDecoda[2,],omnibus[2,],edgeR[2,],DESeq2[2,])

simlist1=simlist

simlist[is.na(simlist)]=0
simlist=simlist*100
rownames(simlist)=c("power.mbDecoda", "power.omnibus", "power.edgeR", "power.DESeq2",
                    "fdr.mbDecoda", "fdr.omnibus", "fdr.edgeR", "fdr.DESeq2")



##########box plot#########
f=function(i){
  index=sim.seed[,i]
  data.frame(class=rep(c(rep("power",n.sim),rep("empirical FDR",n.sim)),each=4),
             method=rep(rep(c("mbDecoda", "omnibus", "Wrench+edgeR", "Wrench+DESeq2"),each=n.sim),2),
             value=as.vector(t(simlist[,index])))
}

power.fdr=foreach(i=1:8, .combine=rbind) %do% f(i)

out=data.frame(zi=rep(c(0,0.1,0.3,0.5),each=16*n.sim),
               q=rep(rep(c("0.5","5"),4),each=8*n.sim),
               power.fdr)
colnames(out)=c("zi.prob","q","class","method","value")
out$zi.prob=factor(out$zi.prob, levels =c("0","0.1","0.3","0.5"))
out$class=factor(out$class, levels =c("power","empirical FDR"))
out$q=factor(out$q, levels =c("0.5","5"))
out$method=factor(out$method, levels =c("mbDecoda", "omnibus", "Wrench+edgeR", "Wrench+DESeq2"))
line=data.frame(class=factor(c("power","empirical FDR"), levels =c("power","empirical FDR")),y=c(NA,5))

p=ggplot(out, aes(x=zi.prob, y=value, fill=method))+ theme_bw()
p=p+geom_boxplot(aes(col=method))+
  facet_grid(vars(class),vars(q),scales="free",space="free", labeller = labeller(.cols = label_both))+
  geom_hline(data= line, aes(yintercept=y),linetype = "dashed")+ theme(legend.position = "top")+
  labs( y = 'empirical FDR and power (%)',x=expression(eta))+ scale_fill_npg()+scale_color_npg()

p
