library(ggplot2)
library("ggsci")
library(doSNOW)
source("mbDecoda.R")
source("ZINB_sim.R")


ncores=70
n.sim=100
K=c(20,50,100,100,100,100,200)
n=c(50,50,30,50,100,200,50)
sig.prob = 0.2
bias = "small"
sim.seed=matrix(1:(n.sim*7),n.sim,7)
confounder = F

DATA=list()
for (i in 1:7) {
  K.i=K[i]
  n1.i=n[i]/2
  n2.i=n[i]/2
  for (j in sim.seed[,i]) {
    set.seed(j)
    DATA[[j]]=ZINB_sim(seed=j,K=K.i,n1=n1.i,n2=n2.i,p=sig.prob,bias =bias,zi=0.3,confounder =confounder)
  }
}


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
simlist=foreach(i = 1:length(DATA), .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(MicrobiomeStat)
  library(glmmTMB)
  library(fastANCOM)
  
  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["confounder"]]
  meta=data.frame(Sample.ID=rownames(count), group)
  OTU = otu_table(t(count), taxa_are_rows = TRUE)
  META = sample_data(meta)
  PHYSEQ = phyloseq(OTU, META)
  
  suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                              p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                      silent = TRUE))
  if (inherits(out_ANCOMBC, "try-error")) {
    fdr_ANCOMBC=0; power_ANCOMBC=0
  }else{
    qval <- out_ANCOMBC[["res"]][["q_val"]]
    id_ANCOMBC <- which(qval <= 0.05)
    power_ANCOMBC=length(intersect(id_ANCOMBC,id))/length(id)
    fdr_ANCOMBC=1-length(intersect(id_ANCOMBC,id))/length(id_ANCOMBC)
  }
  
  suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
  if (inherits(out_ANCOM, "try-error")) {
    fdr_ANCOM=0; power_ANCOM=0
  }else{
    id_ANCOM=which(out_ANCOM[["results"]][["final"]][["REJECT"]]==T)
    power_ANCOM=length(intersect(id_ANCOM,id))/length(id)
    fdr_ANCOM=1-length(intersect(id_ANCOM,id))/length(id_ANCOM)
  }
  
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
  
  #LinDA
  suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
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
cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
locom=foreach(i = 1:length(DATA), .combine = 'cbind') %dopar% {
  library(LOCOM)
  
  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                    , silent = TRUE))
  
  if (inherits(out_LOCOM, "try-error")) {
    fdr_LOCOM=0; power_LOCOM=0
  }else if(is.na(out_LOCOM)){
    fdr_LOCOM=NA; power_LOCOM=0
  }else{
    id_LOCOM=which(out_LOCOM[["q.otu"]]<0.05)
    power_LOCOM=length(intersect(id_LOCOM,id))/length(id)
    fdr_LOCOM=1-length(intersect(id_LOCOM,id))/length(id_LOCOM)
  }
  c(power_LOCOM,fdr_LOCOM)
}
stopCluster(cl)


simlist=rbind(simlist[1:4,],locom[1,],simlist[5:8,],locom[2,])

simlist[is.na(simlist)]=0
simlist=simlist*100
rownames(simlist)=c("power.mbDecoda", "power.ANCOMBC", "power.fastANCOM", "power.LinDA", "power.LOCOM",
                    "fdr.mbDecoda", "fdr.ANCOMBC", "fdr.fastANCOM", "fdr.LinDA", "fdr.LOCOM")

rowMeans(simlist)

