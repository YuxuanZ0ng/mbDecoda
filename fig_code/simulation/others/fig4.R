library(ggplot2)
library("ggsci")
library(stats)
library(MASS)
library(doSNOW)
source("mbDecoda.R")

ncore=40

DATA_ZIP=readRDS(file = "DATA_ZIP_sig0.2.RDS")
DATA_DAS=readRDS(file = "DATA_DAS_sig0.2.RDS")
DATA_LINDA=readRDS(file = "DATA_LINDA_sig0.2.RDS")
DATA_ANCOMBC=readRDS(file = "DATA_ANCOMBC_sig0.2.RDS")





#ZIP
cl <- makeCluster(ncore, type = "SOCK") 
registerDoSNOW(cl)
simlist_ZIP=foreach(i = 1:length(DATA_ZIP), .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(MicrobiomeStat)
  library(glmmTMB)
  library(fastANCOM)
  
  adjust="BH"
  data=DATA_ZIP[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  Gamma=NULL
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
locom=foreach(i = 1:length(DATA_ZIP), .combine = 'cbind') %do% {
  library(LOCOM)
  
  adjust="BH"
  data=DATA_ZIP[[i]]
  id=which(data[["diff.taxa"]]==1)
  group=data[["grp"]]
  count=data[["count"]]
  
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = ncore)
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


simlist_ZIP=rbind(simlist_ZIP[1:4,],locom[1,],simlist_ZIP[5:8,],locom[2,])

simlist_ZIP[is.na(simlist_ZIP)]=0
rownames(simlist_ZIP)=c("power.mbDecoda", "power.ANCOMBC", "power.fastANCOM", "power.LinDA", "power.LOCOM",
                        "fdr.mbDecoda", "fdr.ANCOMBC", "fdr.fastANCOM", "fdr.LinDA", "fdr.LOCOM")




#DAS
cl <- makeCluster(ncore, type = "SOCK") 
registerDoSNOW(cl)
simlist_DAS=foreach(i = DATA_DAS, .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(MicrobiomeStat)
  library(glmmTMB)
  library(fastANCOM)
  
  data=i@otu_table
  count=data@.Data
  adjust="BH"
  otu.tab=t(count)
  id=grep("[[:print:]]+\\-TP$", taxa_names(data))
  grp=factor(gsub("[[:print:]]+\\;", "", sample_names(data)))
  group=as.numeric(grp)-1
  Gamma=NULL
  meta=data.frame(Sample.ID=colnames(otu.tab), group)
  OTU = otu_table(otu.tab, taxa_are_rows = TRUE)
  META = sample_data(meta)
  rownames(META)=colnames(otu.tab)
  PHYSEQ = phyloseq(OTU, META)
  
  #ANCOMBC
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
  
  #fastANCOM
  suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut=1), silent = TRUE))
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
locom=foreach(i = DATA_DAS, .combine = 'cbind') %do% {
  library(LOCOM)
  
  adjust="BH"
  data=i@otu_table
  count=data@.Data
  id=grep("[[:print:]]+\\-TP$", taxa_names(data))
  grp=factor(gsub("[[:print:]]+\\;", "", sample_names(data)))
  group=as.numeric(grp)-1
  
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 50)
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


simlist_DAS=rbind(simlist_DAS[1:4,],locom[1,],simlist_DAS[5:8,],locom[2,])


simlist_DAS[is.na(simlist_DAS)]=0
rownames(simlist_DAS)=c("power.mbDecoda", "power.ANCOMBC", "power.fastANCOM", "power.LinDA", "power.LOCOM",
                        "fdr.mbDecoda", "fdr.ANCOMBC", "fdr.fastANCOM", "fdr.LinDA", "fdr.LOCOM")
rowMeans(simlist_DAS)



#LINDA
cl <- makeCluster(ncore, type = "SOCK") 
registerDoSNOW(cl)
simlist_LINDA=foreach(i = DATA_LINDA, .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(MicrobiomeStat)
  library(glmmTMB)
  library(fastANCOM)
  
  adjust="BH"
  data=i
  id=which(data[["H"]]==1)
  group=as.numeric(data[["Z"]][["u"]])-1
  otu.tab=data[["Y"]]
  meta=data.frame(Sample.ID=colnames(otu.tab), group)
  OTU = otu_table(otu.tab, taxa_are_rows = TRUE)
  META = sample_data(meta)
  rownames(META)=colnames(otu.tab)
  PHYSEQ = phyloseq(OTU, META)
  count <-t(otu.tab)
  
  #ANCOMBC
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
  
  #fastANCOM
  suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut=1), silent = TRUE))
  if (inherits(out_ANCOM, "try-error")) {
    fdr_ANCOM=0; power_ANCOM=0
  }else{
    id_ANCOM=which(out_ANCOM[["results"]][["final"]][["REJECT"]]==T)
    power_ANCOM=length(intersect(id_ANCOM,id))/length(id)
    fdr_ANCOM=1-length(intersect(id_ANCOM,id))/length(id_ANCOM)
  }
  
  
  #mbDecoda
  suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=NULL, W=NULL, adjust=adjust)
                                    , silent = TRUE))
  if (inherits(out_mbDecoda, "try-error")) {
    fdr_mbDecoda=0; power_mbDecoda=0
  }else{
    id_mbDecoda=which(out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]<0.05)
    power_mbDecoda=length(intersect(id_mbDecoda,id))/length(id)
    fdr_mbDecoda=1-length(intersect(id_mbDecoda,id))/length(id_mbDecoda)
  }
  
  suppressWarnings(out_LinDA <- try(linda(otu.tab, meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
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
locom=foreach(i = DATA_LINDA, .combine = 'cbind') %do% {
  library(LOCOM)
  
  adjust="BH"
  data=i
  id=which(data[["H"]]==1)
  group=as.numeric(data[["Z"]][["u"]])-1
  otu.tab=data[["Y"]]
  count <-t(otu.tab)
  
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 50)
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


simlist_LINDA=rbind(simlist_LINDA[1:4,],locom[1,],simlist_LINDA[5:8,],locom[2,])


simlist_LINDA[is.na(simlist_LINDA)]=0
rownames(simlist_LINDA)=c("power.mbDecoda", "power.ANCOMBC", "power.fastANCOM", "power.LinDA", "power.LOCOM",
                          "fdr.mbDecoda", "fdr.ANCOMBC", "fdr.fastANCOM", "fdr.LinDA", "fdr.LOCOM")
rowMeans(simlist_LINDA)



#ANCOMBC
cl <- makeCluster(ncore, type = "SOCK") 
registerDoSNOW(cl)

simlist_ANCOMBC=foreach(i = DATA_ANCOMBC, .combine = 'cbind') %dopar% {
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(MicrobiomeStat)
  library(glmmTMB)
  library(fastANCOM)
  
  test.dat=i
  obs.abn=test.dat$obs.abn
  # Pre-processing
  feature.table=obs.abn
  adjust="BH"
  group=test.dat[["grp"]]-1
  count <-t(feature.table)
  meta=data.frame(Sample.ID=colnames(feature.table), group)
  OTU = otu_table(feature.table, taxa_are_rows = TRUE)
  META = sample_data(meta)
  rownames(META)=colnames(feature.table)
  PHYSEQ = phyloseq(OTU, META)
  id=which(test.dat[["effect.size"]]!=1)
  
  #ANCOMBC
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
  
  #fastANCOM
  suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
  if (inherits(out_ANCOM, "try-error")) {
    fdr_ANCOM=0; power_ANCOM=0
  }else{
    id_ANCOM=which(out_ANCOM[["results"]][["final"]][["REJECT"]]==T)
    power_ANCOM=length(intersect(id_ANCOM,id))/length(id)
    fdr_ANCOM=1-length(intersect(id_ANCOM,id))/length(id_ANCOM)
  }
  
  
  #mbDecoda
  suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=NULL, W=NULL, adjust=adjust)
                                    , silent = TRUE))
  if (inherits(out_mbDecoda, "try-error")) {
    fdr_mbDecoda=0; power_mbDecoda=0
  }else{
    id_mbDecoda=which(out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]<0.05)
    power_mbDecoda=length(intersect(id_mbDecoda,id))/length(id)
    fdr_mbDecoda=1-length(intersect(id_mbDecoda,id))/length(id_mbDecoda)
  }
  
  suppressWarnings(out_LinDA <- try(linda(feature.table, meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
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
locom=foreach(i = DATA_ANCOMBC, .combine = 'cbind') %do% {
  library(LOCOM)
  
  test.dat=i
  obs.abn=test.dat$obs.abn
  adjust="BH"
  group=test.dat[["grp"]]-1
  count <-t(obs.abn)
  id=which(test.dat[["effect.size"]]!=1)
  
  suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 50)
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


simlist_ANCOMBC=rbind(simlist_ANCOMBC[1:4,],locom[1,],simlist_ANCOMBC[5:8,],locom[2,])


simlist_ANCOMBC[is.na(simlist_ANCOMBC)]=0
rownames(simlist_ANCOMBC)=c("power.mbDecoda", "power.ANCOMBC", "power.fastANCOM", "power.LinDA", "power.LOCOM",
                            "fdr.mbDecoda", "fdr.ANCOMBC", "fdr.fastANCOM", "fdr.LinDA", "fdr.LOCOM")
rowMeans(simlist_ANCOMBC)



simlist=rbind(simlist_ZIP,simlist_DAS,simlist_LINDA,simlist_ANCOMBC)
out=data.frame(simulation=rep(c("M1","M4","M3","M2"),each=10),
               class=rep(rep(c("power","empirical FDR"),each=5),4),
               method=rep(c("mbDecoda","ANCOMBC","fastANCOM","LinDA","LOCOM"),8),
               value=rowMeans(simlist*100))
out$simulation=factor(out$simulation, levels =c("M1","M2","M3","M4"))
out$class=factor(out$class, levels =c("power","empirical FDR"))
out$method=factor(out$method, levels =c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda"))
out$value=as.numeric(out$value)

line=data.frame(class=factor(c("power","empirical FDR"), levels =c("power","empirical FDR")),y=c(NA,5))

p<-ggplot(out,aes(method,value,color=method,fill=method))+
  scale_y_continuous(expand = c(0,0))+
  geom_bar(stat="identity",position=position_dodge(0.75),width = 0.5)+ #??״ͼ
  facet_grid(class~simulation,scales = 'free',space = "free_y")+#????
  geom_hline(data= line, aes(yintercept=y),linetype = "dashed")+
  scale_x_discrete(guide ="none")+
  labs( y = 'empirical FDR and power (%)')+theme_bw()+ 
  scale_fill_npg()+scale_color_npg()+ theme(legend.position = "top")
p

