library(dplyr)
library(mia)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(phyloseq)
library(doSNOW)
library(curatedMetagenomicData)
library(MASS)
source("mbDecoda.R")


nsim=100
ncores=50

data.gen=function(count, Gamma){
  group=rbinom(length(count[,1]),1,0.5)
  test.data=list(count, group, Gamma)
  names(test.data)=c("count", "grp", "Gamma")
  return(test.data)
}


################## DATA #################
HMP_2012=curatedMetagenomicData("HMP_2012.relative_abundance", dryrun = FALSE, counts = TRUE)%>%
  mergeData()
altExps(HMP_2012) <-
  splitByRanks(HMP_2012)
HMP_2012_class=altExp(HMP_2012,"class")
otutable=HMP_2012_class@assays@data@listData[["relative_abundance"]]
otutable=t(as(otutable,"matrix"))


datasize <- nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  0.05)
otutable = otutable[,keepOTUs]

as.numeric(colSums(otutable==0)/748)


nsim=100
ncores=50
Gamma=NULL

DATA=list()
for (i in 1:nsim) {
  set.seed(i)
  DATA[[i]]=data.gen(count=otutable,Gamma=Gamma)
}




################## simulation #################
cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_ANCOMBC=foreach(i = DATA, .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  
  adjust="BH"
  data=i
  group=data[["grp"]]
  count=data[["count"]]
  meta=data.frame(Sample.ID=rownames(count), group)
  rownames(meta)=rownames(count)
  OTU = otu_table(t(count), taxa_are_rows = TRUE)
  META = sample_data(meta)
  PHYSEQ = phyloseq(OTU, META)
  
  
  suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                              p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                      silent = TRUE))
  if (inherits(out_ANCOMBC, "try-error")) {
    p=NULL
    p1=NULL
  }else{
    p <- out_ANCOMBC[["res"]][["q_val"]][["group"]]
    p1 <- out_ANCOMBC[["res"]][["p_val"]][["group"]]
  }
  c(p,p1)
}
stopCluster(cl)


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_ANCOM=foreach(i = DATA, .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(fastANCOM)
  
  adjust="BH"
  data=i
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["Gamma"]]
  
  
  suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, Z=Gamma, zero_cut =1), silent = TRUE))
  if (inherits(out_ANCOM, "try-error")) {
    p = NULL
    p1 = NULL
  }else{
    p = out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
    p1 = out_ANCOM[["results"]][["final"]][["log2FC.pval"]]
  }
  c(p,p1)
}
stopCluster(cl)


cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_mbDecoda=foreach(i = DATA, .combine = 'cbind') %dopar% {
  library(glmmTMB)
  
  adjust="BH"
  data=i
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["Gamma"]]
  
  
  #mbDecoda
  suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, interval.prob=0.9, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                       , silent = TRUE))
  if (inherits(out_mbDecoda, "try-error")) {
    p = NULL
    p1 = NULL
  }else{
    p = out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
    p1 = out_mbDecoda[["Bias_correct"]][["DAA"]][["p.val"]]
  }
  c(p,p1)
  
}
stopCluster(cl)



cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_LinDA=foreach(i = DATA, .combine = 'cbind') %dopar% {
  library(MicrobiomeStat)
  
  adjust="BH"
  data=i
  group=data[["grp"]]
  count=data[["count"]]
  meta=data.frame(Sample.ID=rownames(count), group)
  rownames(meta)=rownames(count)
  
  #LinDA
  suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                    , silent = TRUE))
  if (inherits(out_LinDA, "try-error")) {
    p = NULL
    p1 = NULL
  }else{
    p = out_LinDA[["output"]][["group"]][["padj"]]
    p1 = out_LinDA[["output"]][["group"]][["pvalue"]]
  }
  c(p,p1)
  
}
stopCluster(cl)

cl <- makeCluster(ncores, type = "SOCK") 
registerDoSNOW(cl)
p_LOCOM=foreach(i = DATA, .combine = 'cbind') %dopar% {
  library(LOCOM)
  
  adjust="BH"
  data=i
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=data[["Gamma"]]
  
  
  #LOCOM
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), C=Gamma, fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                    , silent = TRUE))
  if (inherits(out_LOCOM, "try-error")) {
    p = NULL
    p1 = NULL
  }else{
    p = t(out_LOCOM[["q.otu"]])
    p1 = t(out_LOCOM[["p.otu"]])
  }
  c(p,p1)
}
stopCluster(cl)



HMP_class=data.frame(rowSums(p_ANCOMBC<0.05)/100,rowSums(p_ANCOM<0.05)/100,rowSums(p_LinDA<0.05)/100,rowSums(p_LOCOM<0.05)/100,rowSums(p_mbDecoda<0.05)/100)
HMP_class[is.na(HMP_class)]=0
colnames(HMP_class)=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
HMP_class_q=HMP_class[1:18,]
HMP_class_q$taxon=1:18
rownames(HMP_class_q)=colnames(DATA2[[1]][["count"]])
HMP_class_p=HMP_class[19:36,]
HMP_class_p$taxon=1:18
rownames(HMP_class_p)=colnames(DATA2[[1]][["count"]])


HMP_class_long <- pivot_longer(HMP_class_q, cols = c("ANCOMBC", "fastANCOM", "LinDA", "LOCOM", "mbDecoda"), 
                               names_to = "Method", values_to = "Value")

# 绘制分面散点图
ggplot(HMP_class_long, aes(x = factor(taxon), y = Value)) + 
  geom_point() + 
  facet_wrap(~ Method, scales = "fixed") + ylim(c(0,0.2))+
  labs(x = "Taxon", y = "type I error") + theme_bw()+ 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+ 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") 


HMP_class_long1 <- pivot_longer(HMP_class_p, cols = c("ANCOMBC", "fastANCOM", "LinDA", "LOCOM", "mbDecoda"), 
                                names_to = "Method", values_to = "Value")

# 绘制分面散点图
ggplot(HMP_class_long1, aes(x = factor(taxon), y = Value)) + 
  geom_point() + 
  facet_wrap(~ Method, scales = "fixed") + ylim(c(0,0.2))+
  labs(x = "Taxon", y = "type I error") + theme_bw()+ 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+ 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") 
