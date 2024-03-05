library(curatedMetagenomicData)
library(dplyr)
library(mia)
library(tidyverse)
library(phyloseq)
library(doSNOW)
library(MASS)
library(truncnorm)
source("mbDecoda.R")


###########################case 1############################
ncores=50
n.sim=100
K=100
n1=25
n2=25
sig.prob = 0.2
bias = "small"
sim.seed=1:n.sim
confounder = F


ZINB_sim1=function(seed, K, n1, n2, p, bias=c("large", "small"), overdisp=3, confounder=c("TRUE", "FALSE")){
  # Total number of samples
  n=n1+n2
  x=c(rep(0,n1),rep(1,n2))
  w=c(0.6,0.3,0.1)
  
  # generate baseline beta0，ebeta0 means exp(beta0)
  set.seed(seed)
  index=sample(c(1, 2, 3), K, replace = T, prob=w) 
  ebeta0=rep(NA, K)
  set.seed(seed)
  ebeta0[which(index==1)]=rgamma(length(which(index==1)), shape=50, rate=1)
  ebeta0[which(index==2)]=rgamma(length(which(index==2)), shape=200, rate=1)
  ebeta0[which(index==3)]=rgamma(length(which(index==3)), shape=10000, rate=1)
  beta0=log(ebeta0)
  
  # generate H
  set.seed(seed)
  H=rep(0,K)
  C1=1:floor(K*p)
  wt=runif(1, 0, 1)
  C2=sample(C1, wt*length(C1), replace=FALSE)
  H[C1]=1
  
  # generate delta, edelta means exp(delta)
  set.seed(seed)
  edelta=rep(1, K)
  edelta[C1]=runif(length(C1), 1.5, 10)
  edelta[C2]=1/runif(length(C2), 1.5, 10)
  delta=log(edelta)
  
  # generate d
  if(bias=="small"){
    d_bias=mean(abs(delta[C1]))/2
  }else{
    d_bias=mean(abs(delta[C1]))*5
  }
  set.seed(seed)
  d=c(rnorm(n1,0,1),rnorm(n2,d_bias,1))
  
  # generate Y, lY means log(Y)
  lY=matrix(NA,n,K)
  for (i in 1:n) {
    lY[i,]=d[i]+beta0+x[i]*delta
  }
  #confounders
  if(confounder){
    set.seed(seed)
    Gamma=cbind(runif(n,0,1),rbinom(n,1,0.5))
    beta=rbind(rep(-1,K),rep(1,K))
    lY=lY+Gamma%*%beta
    rownames(Gamma)=paste0("sub", seq(n))
  }else{Gamma=NULL}
  Y=exp(lY)
  
  #generate Z
  set.seed(seed)
  eta <- runif(K,0,0.7)
  C3=sample(1:K, 0.5*K, replace=FALSE)
  eta[C3]=0
  Z <- matrix(0,n,K)
  set.seed(seed)
  for(j in 1:n){
    Z[j,] <- rbinom(K, size=1, prob=eta)
  }
  
  #generate C
  C=matrix(NA, ncol=K, nrow=n)
  set.seed(seed)
  for(i in 1:n){
    for(k in 1:K){
      C[i,k]=rnegbin(1,Y[i,k],theta =overdisp)
    }
  }
  C[Z==1]=0
  C=as.data.frame(C)
  Z=as.data.frame(Z)
  
  # Prepare outputs
  rownames(Z)=paste0("sub", seq(n))
  colnames(Z)=paste0("taxon", seq(K))
  rownames(C)=paste0("sub", seq(n))
  colnames(C)=paste0("taxon", seq(K))
  names(H)=paste0("taxon", seq(K))
  names(x)=paste0("sub", seq(n))
  names(delta)=paste0("taxon", seq(K))
  names(d)=paste0("sub", seq(n))
  
  test.data=list(C, delta, x, Gamma, H, d, Z, eta)
  names(test.data)=c("count", "delta", "grp", "confounder", "diff.taxa", "d", "zero.pos","zero.prop")
  return(test.data)
}

DATA=list()
for (j in sim.seed) {
  set.seed(j)
  DATA[[j]]=ZINB_sim1(seed=j, K=K, n1=n1, n2=n2, p=sig.prob, bias =bias, confounder=confounder)
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
locom=foreach(i = 1:length(DATA), .combine = 'cbind') %do% {
  library(LOCOM)
  
  adjust="BH"
  data=DATA[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 40)
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


simlist=rbind(simlist[1:4,],locom[1,],simlist[5:8,],locom[2,])

simlist[is.na(simlist)]=0
simlist=simlist*100
rownames(simlist)=c("power.mbDecoda", "power.ANCOMBC", "power.fastANCOM", "power.LinDA", "power.LOCOM",
                    "fdr.mbDecoda", "fdr.ANCOMBC", "fdr.fastANCOM", "fdr.LinDA", "fdr.LOCOM")
rowMeans(simlist)

