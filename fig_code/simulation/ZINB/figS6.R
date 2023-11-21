library(ggplot2)
library("ggsci")
source("mbDecoda.R")
source("ZINB_sim.R")

library(doSNOW)
n.sim=100
K=rep(100,4)
n1=rep(25,4)
n2=n1
sig.prob=0.2
zi.prob = c(0,0.1,0.3,0.5)
bias = "small"
sim.seed=matrix(1:(n.sim*4),n.sim,4)
confounder = T
ncores=40


DATA=list()
for (i in 1:4) {
  zi.prob.i=zi.prob[i]
  K.i=K[i]
  n1.i=n1[i]
  n2.i=n2[i]
  for (j in sim.seed[,i]) {
    set.seed(j)
    DATA[[j]]=ZINB_sim(seed=j,K=K.i,n1=n1.i,n2=n2.i,p=sig.prob,bias =bias,zi=zi.prob.i,confounder =confounder)
  }
}




cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
simlist=foreach(i = DATA, .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(MicrobiomeStat)
  library(glmmTMB)
  library(fastANCOM)
  
  adjust="BH"
  data=i
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["confounder"]]
  meta=data.frame(Sample.ID=rownames(count), group,x1=Gamma[,1],x2=Gamma[,2])
  OTU = otu_table(t(count), taxa_are_rows = TRUE)
  META = sample_data(meta)
  PHYSEQ = phyloseq(OTU, META)
  zi=data[["zi.prop"]]
  K=length(count[1,])
  n=length(count[,1])
  
  suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group+x1+x2",
                                              p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                      silent = TRUE))
  if (inherits(out_ANCOMBC, "try-error")) {
    fdr_ANCOMBC=0; power_ANCOMBC=0
  }else{
    qval <- out_ANCOMBC[["res"]][["q_val"]][["group"]]
    id_ANCOMBC <- which(qval <= 0.05)
    power_ANCOMBC=length(intersect(id_ANCOMBC,id))/length(id)
    fdr_ANCOMBC=1-length(intersect(id_ANCOMBC,id))/length(id_ANCOMBC)
  }
  
  suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, Z=Gamma, zero_cut =1), silent = TRUE))
  if (inherits(out_ANCOM, "try-error")) {
    fdr_ANCOM=0; power_ANCOM=0
  }else{
    id_ANCOM=which(out_ANCOM[["results"]][["final"]][["REJECT"]]==T)
    power_ANCOM=length(intersect(id_ANCOM,id))/length(id)
    fdr_ANCOM=1-length(intersect(id_ANCOM,id))/length(id_ANCOM)
  }
  
  #mbDecoda
  suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, adjust=adjust)
                                    , silent = TRUE))
  if (inherits(out_mbDecoda, "try-error")) {
    fdr_mbDecoda=0; power_mbDecoda=0
  }else{
    id_mbDecoda=which(out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]<0.05)
    power_mbDecoda=length(intersect(id_mbDecoda,id))/length(id)
    fdr_mbDecoda=1-length(intersect(id_mbDecoda,id))/length(id_mbDecoda)
  }
  
  #LinDA
  suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group+x1+x2',alpha = 0.05,p.adj.method = adjust)
                                    , silent = TRUE))
  if (inherits(out_LinDA, "try-error")) {
    fdr_LinDA=0; power_LinDA=0
  }else{
    id_LinDA=which(out_LinDA[["output"]][["group"]][["padj"]]<0.05)
    power_LinDA=length(intersect(id_LinDA,id))/length(id)
    fdr_LinDA=1-length(intersect(id_LinDA,id))/length(id_LinDA)
  }
  

  
  c(power_mbDecoda, power_ANCOMBC, power_ANCOM, power_LinDA,
    fdr_mbDecoda, fdr_ANCOMBC, fdr_ANCOM, fdr_LinDA)
}

stopCluster(cl)



# LOCOM
locom=foreach(i = 1:length(DATA), .combine = 'cbind') %do% {
  library(LOCOM)

  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["confounder"]]
  
  # LOCOM
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), C=Gamma, fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = ncores)
                                    , silent = TRUE))
  if (inherits(out_LOCOM, "try-error")) {
    fdr_LOCOM=0; power_LOCOM=0
  }else{
    id_LOCOM=which(out_LOCOM[["q.otu"]]<0.05)
    power_LOCOM=length(intersect(id_LOCOM,id))/length(id)
    fdr_LOCOM=1-length(intersect(id_LOCOM,id))/length(id_LOCOM)
  }
  
  c(power_LOCOM,fdr_LOCOM)
}


simlist=rbind(simlist[1:4,],locom[1,],simlist[5:8,],locom[2,])

simlist[is.na(simlist)]=0
simlist=simlist*100
rownames(simlist)=c("power.mbDecoda", "power.ANCOMBC", "power.fastANCOM", "power.LinDA", "power.LOCOM", "fdr.mbDecoda", "fdr.ANCOMBC", "fdr.fastANCOM", "fdr.LinDA", "fdr.LOCOM")


f=function(i){
  index=sim.seed[,i]
  data.frame(class=rep(c("power","empirical FDR"),each=5),
             method=rep(c("mbDecoda", "ANCOMBC", "fastANCOM", "LinDA", "LOCOM"),2),
             value=rowMeans(simlist[,index]))
}

power.fdr=foreach(i=1:4, .combine=rbind) %do% f(i)

out=data.frame(zi=rep(c(0,0.1,0.3,0.5),each=10),
               power.fdr)
colnames(out)=c("zi.prob","class","method","value")
out$zi.prob=factor(out$zi.prob, levels =c("0","0.1","0.3","0.5"))
out$class=factor(out$class, levels =c("power","empirical FDR"))
line=data.frame(class=factor(c("power","empirical FDR"), levels =c("power","empirical FDR")),y=c(NA,5))

p=ggplot(out, aes(x=zi.prob, y=value, fill=method))+ theme_bw()
p=p+geom_bar(aes(col=method),stat="identity",position=position_dodge(0.75),width = 0.5)+
  facet_grid(vars(class), labeller = labeller(.cols = label_both))+
  geom_hline(data= line, aes(yintercept=y),linetype = "dashed")+ theme(legend.position = "top")+
  labs( y = 'empirical FDR and power (%)',x=expression(eta))+ scale_fill_npg()+scale_color_npg()

p

