library(ggplot2)
library(ggsci)
library(doSNOW)
library(tidyverse)
library(phyloseq)
library(ANCOMBC)
library(MicrobiomeStat)
library(glmmTMB)
library(fastANCOM)
library(LOCOM)
library(UpSetR)
library(ggstar)
library(aplot)
library(ggVennDiagram)
library(patchwork)
library(vegan)
library(RColorBrewer)
library(scales)
library(ggthemes)
library(dplyr)
library(ggrepel)


source("mbDecoda.R")

CRC_baxter=readRDS("CRC_baxter.RDS")
CRC_zackular=readRDS("CRC_zackular.RDS")

CDI_schubert=readRDS("CDI_schubert.RDS")
CDI_vincent=readRDS("CDI_vincent.RDS")

ME_NagySzakalD=readRDS("ME_NagySzakalD.RDS")
ME_NagySzakalD_sub=readRDS("ME_NagySzakalD_sub.RDS")

IGT_KarlssonFH=readRDS("IGT_KarlssonFH.RDS")
IGT_KarlssonFH_sub=readRDS("IGT_KarlssonFH_sub.RDS")



##########################NMDS#######################
#####CRC#####
otu=rbind(CRC_baxter[[1]],CRC_zackular[[1]])
otu=otu/rowSums(otu)
otu %>%data.frame() ->otu
group=c(CRC_baxter[[2]],CRC_zackular[[2]]+2)
group[group==0]="Baxter:H"
group[group==1]="Baxter:CRC"
group[group==2]="Zackular:H"
group[group==3]="Zackular:CRC"
distance <- vegdist(otu, method = 'bray')
dfNmds<-metaMDS(otu,distance="bray",k = 2)
plot_data <- data.frame(dfNmds$points)
plot_data$group=as.factor(group)
names(plot_data)[1:2] <- c('NMDS1', 'NMDS2')
centroid1 <- aggregate(cbind(NMDS1,NMDS2) ~ group,
                      data = plot_data,
                      FUN = mean)
centroid1$reference=rep(c("Baxter","Zackular"),each=2)
centroid1$X=rep(c("case","control"),2)
dist(centroid1[,2:3],p=2)
p1=ggplot(data = centroid1, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(x=NMDS1, y=NMDS2,color=as.factor(reference)),size=5) +
  geom_text_repel(aes(label = group))+
  xlim(c(-0.08,0.18))+ylim(c(-0.08,0.18))+
  geom_line(aes(group = X),color= "#5A9599FF")+
  geom_line(aes(group = reference),color= c("#FF6348FF"))+
  scale_color_manual(values = c("#D55E0080","#91D1C2FF"))+
  labs(title = expression("CRC"))+
  theme_bw()+theme(plot.title = element_text(hjust=0.5),legend.position = "none")
p1

#####CDI#####
otu=rbind(CDI_schubert[[1]],CDI_vincent[[1]])
otu=otu/rowSums(otu)
otu %>%data.frame() ->otu
group=c(CDI_schubert[[2]],CDI_vincent[[2]]+2)
group[group==0]="Schubert:H"
group[group==1]="Schubert:CDI"
group[group==2]="Vincent:H"
group[group==3]="Vincent:CDI"
distance <- vegdist(otu, method = 'bray')
dfNmds<-metaMDS(otu,distance="bray",k = 2)
plot_data <- data.frame(dfNmds$points)
plot_data$group=as.factor(group)
names(plot_data)[1:2] <- c('NMDS1', 'NMDS2')
centroid2 <- aggregate(cbind(NMDS1,NMDS2) ~ group,
                      data = plot_data,
                      FUN = mean)

centroid2$reference=rep(c("Schubert","Vincent"),each=2)
centroid2$X=rep(c("case","control"),2)
p2=ggplot(data = centroid2, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(x=NMDS1, y=NMDS2,color=as.factor(reference)),size=5) +
  geom_text_repel(aes(label = group))+
  xlim(c(-0.9,0.7))+ylim(c(-0.9,0.7))+
  geom_line(aes(group = X),color= "#5A9599FF")+
  geom_line(aes(group = reference),color= c("#FF6348FF"))+
  scale_color_manual(values = c("#D55E0080","#91D1C2FF"))+
  labs(title = expression("CDI"))+
  theme_bw()+theme(plot.title = element_text(hjust=0.5),legend.position = "none")
p2


####IGT#####
otu=rbind(IGT_KarlssonFH[[1]],IGT_KarlssonFH_sub[[1]])
otu=otu/rowSums(otu)
otu %>%data.frame() ->otu
group=c(IGT_KarlssonFH[[2]],IGT_KarlssonFH_sub[[2]]+2)
group[group==0]="Full:H"
group[group==1]="Full:IGT"
group[group==2]="Sub:H"
group[group==3]="Sub:IGT"
group=as.factor(group)
distance <- vegdist(otu, method = 'bray')
dfNmds<-metaMDS(otu,distance="bray",k = 2)
plot_data <- data.frame(dfNmds$points)
plot_data$group=as.factor(group)
names(plot_data)[1:2] <- c('NMDS1', 'NMDS2')
centroid3 <- aggregate(cbind(NMDS1,NMDS2) ~ group,
                      data = plot_data,
                      FUN = mean)

centroid3$scale=rep(c("Full","Sub"),each=2)
centroid3$X=rep(c("case","control"),2)
p3=ggplot(data = centroid3, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(x=NMDS1, y=NMDS2,color=scale),size=5) +
  geom_text_repel(aes(label = group))+
  xlim(c(-0.15,0.15))+ylim(c(-0.15,0.15))+
  geom_line(aes(group = X),color= "#5A9599FF")+
  geom_line(aes(group = scale),color= c("#FF6348FF"))+
  scale_color_manual(values = c("#D55E0080","#91D1C2FF"))+
  labs(title = expression("IGT"))+
  theme_bw()+theme(plot.title = element_text(hjust=0.5),legend.position = "none")
p3




####ME/CFS#####
otu=rbind(ME_NagySzakalD[[1]],ME_NagySzakalD_sub[[1]])
otu=otu/rowSums(otu)
otu %>%data.frame() ->otu
group=c(ME_NagySzakalD[[2]],ME_NagySzakalD_sub[[2]]+2)
group[group==0]="Full:H"
group[group==1]="Full:ME/CFS"
group[group==2]="Sub:H"
group[group==3]="Sub:ME/CFS"
distance <- vegdist(otu, method = 'bray')
dfNmds<-metaMDS(otu,distance="bray",k = 2)
plot_data <- data.frame(dfNmds$points)
plot_data$group=as.factor(group)
names(plot_data)[1:2] <- c('NMDS1', 'NMDS2')
centroid4 <- aggregate(cbind(NMDS1,NMDS2) ~ group,
                      data = plot_data,
                      FUN = mean)

centroid4$scale=rep(c("Full","Sub"),each=2)
centroid4$X=rep(c("case","control"),2)
p4=ggplot(data = centroid4, aes(x=NMDS1, y=NMDS2))+
  geom_point(aes(x=NMDS1, y=NMDS2,color=scale),size=5)+
  geom_text_repel(aes(label = group))+
  xlim(c(-0.15,0.17))+ylim(c(-0.15,0.17))+
  geom_line(aes(group = X),color= "#5A9599FF")+
  geom_line(aes(group = scale),color= c("#FF6348FF"))+
  scale_color_manual(values = c("#D55E0080","#91D1C2FF"))+
  labs(title = expression("ME/CFS"))+
  theme_bw()+theme(plot.title = element_text(hjust=0.5),legend.position = "none")
p4






######################consistance####################
###CRC###
###CRC_baxter###
count=CRC_baxter[[1]]
group=CRC_baxter[[2]]
adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC=NULL
}else{
  q_ANCOMBC <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM=NULL
}else{
  q_ANCOM=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda=NULL
}else{
  q_mbDecoda=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA=NULL
}else{
  q_LinDA=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM=NULL
}else{
  q_LOCOM=t(out_LOCOM[["q.otu"]])
}



###CRC_zackular###
count=CRC_zackular[[1]]
group=CRC_zackular[[2]]
adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC1=NULL
}else{
  q_ANCOMBC1 <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM1=NULL
}else{
  q_ANCOM1=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda1=NULL
}else{
  q_mbDecoda1=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA1=NULL
}else{
  q_LinDA1=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM1=NULL
}else{
  q_LOCOM1=t(out_LOCOM[["q.otu"]])
}

q_sum_CRC=data.frame(q_ANCOMBC,q_ANCOM,q_LinDA,q_LOCOM,q_mbDecoda,q_ANCOMBC1,q_ANCOM1,q_LinDA1,q_LOCOM1,q_mbDecoda1)
k=length(q_sum_CRC)/2


alpha=0.05
f=function(i){
  inter=length(intersect(which(q_sum_CRC[,i]<alpha),which(q_sum_CRC[,i+k]<alpha)))
  consist=inter/length(union(which(q_sum_CRC[,i]<alpha),which(q_sum_CRC[,i+k]<alpha)))
  consist
}


consist=foreach(i=1:k, .combine=rbind) %do% f(i)
consist[is.na(consist)]=0
method=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
out_CRC=data.frame(method,consist)



p5=ggplot(data=out_CRC,aes(method,consist,color=method,fill=method))+
  geom_bar(stat="identity",position=position_dodge(0.75),width = 0.5, color = "transparent")+ ylim(0,0.4)+
  theme_bw()+scale_fill_npg()+
  theme(plot.title = element_text(hjust=0.5),legend.position = "none" , axis.text.x = element_blank())+
  labs(x="method",y="Jaccard",title = expression("CRC"))
p5


findings_0.05=apply((q_sum_CRC<0.05),2,as.numeric)
findings_0.05=as.data.frame(findings_0.05)
rownames(findings_0.05)=rownames(q_sum_CRC)

f1=function(i){
  t=0.05
  inter=rep(0,length(q_sum_CRC[,i]))
  inter[intersect(which(q_sum_CRC[,i]<t),which(q_sum_CRC[,i+k]<t))]=1
  inter
}

consist_0.05=foreach(i=1:k, .combine=cbind) %do% f1(i)
consist_0.05=as.data.frame(consist_0.05)
colnames(consist_0.05)=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
rownames(consist_0.05)=rownames(q_sum_CRC)

df=findings_0.05[rowSums(consist_0.05)==1,]

df=data.frame(paste("taxon", which(rowSums(consist_0.05)==1), sep = ""),df)
colnames(df)=c("OTU","ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda","ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
df1=df[,c(1,1:k+1)]
df2=df[,-(1:k+1)]
df1<-reshape2::melt(df1)
df1$value[df1$value==0]="undetected"
df1$value[df1$value==1]="detected"
df2<-reshape2::melt(df2)
df2$value[df2$value==0]="undetected"
df2$value[df2$value==1]="detected"
df2$value=as.factor(df1$value,levels = c("undetected","detected"))
p6=ggplot()+
  geom_star(data=df1,aes(x=variable,y=OTU,fill=value),size=4,
            starshape=15)+
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x.top = element_text(angle = 90,
                                       hjust = 0,
                                       vjust= 0.5),
        plot.caption = element_text(hjust=0.5),
        legend.position = "none")+
  theme(panel.grid = element_blank())+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "left")+
  labs(caption =expression("Baxter et al. \n n=292"))+
  labs(x=NULL,y=NULL)+
  scale_fill_manual(values = c("#2a83a2","#ffe8d560"))

p7=ggplot()+
  geom_star(data=df2,aes(x=variable,y=OTU,fill=value),size=4,
            starshape=15)+
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x.top = element_text(angle = 90,
                                       hjust = 0,
                                       vjust= 0.5),
        plot.caption = element_text(hjust=0.5),
        legend.position = "right")+
  theme(panel.grid = element_blank())+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "left")+
  labs(caption = expression("Zackular et al. \n n=60"))+
  labs(x=NULL,y=NULL,fill="Status")+
  scale_fill_manual(values = c("#2a83a2","#ffe8d560"))

p8=p6%>%insert_right(p7,width = 1)
p8




###CDI###
###CDI_schubert###
count=CDI_schubert[[1]]
group=CDI_schubert[[2]]
adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC=NULL
}else{
  q_ANCOMBC <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM=NULL
}else{
  q_ANCOM=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda=NULL
}else{
  q_mbDecoda=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA=NULL
}else{
  q_LinDA=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM=NULL
}else{
  q_LOCOM=t(out_LOCOM[["q.otu"]])
}



###CDI_vincent###
count=CDI_vincent[[1]]
group=CDI_vincent[[2]]
adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC1=NULL
}else{
  q_ANCOMBC1 <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM1=NULL
}else{
  q_ANCOM1=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda1=NULL
}else{
  q_mbDecoda1=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA1=NULL
}else{
  q_LinDA1=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM1=NULL
}else{
  q_LOCOM1=t(out_LOCOM[["q.otu"]])
}

q_sum_CDI=data.frame(q_ANCOMBC,q_ANCOM,q_LinDA,q_LOCOM,q_mbDecoda,q_ANCOMBC1,q_ANCOM1,q_LinDA1,q_LOCOM1,q_mbDecoda1)
k=length(q_sum_CDI)/2



f=function(i){
  consist=rep(NA,length(alpha))
  for (j in 1:length(alpha)) {
    t=alpha[j]
    inter=length(intersect(which(q_sum_CDI[,i]<t),which(q_sum_CDI[,i+k]<t)))
    consist[j]=inter/length(union(which(q_sum_CDI[,i]<t),which(q_sum_CDI[,i+k]<t)))
  }
  cbind(alpha,consist)
}


alpha=0.05
f=function(i){
  inter=length(intersect(which(q_sum_CDI[,i]<alpha),which(q_sum_CDI[,i+k]<alpha)))
  consist=inter/length(union(which(q_sum_CDI[,i]<alpha),which(q_sum_CDI[,i+k]<alpha)))
  consist
}


consist=foreach(i=1:k, .combine=rbind) %do% f(i)
consist[is.na(consist)]=0
method=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
out_CDI=data.frame(method,consist)


p9=ggplot(data=out_CDI,aes(method,consist,color=method,fill=method))+
  geom_bar(stat="identity",position=position_dodge(0.75),width = 0.5, color = "transparent")+ ylim(0,0.4)+
  theme_bw()+scale_fill_npg()+
  theme(plot.title = element_text(hjust=0.5),legend.position = "none",axis.text.x = element_blank())+
  labs(x="method",y="Jaccard",title = expression("CDI"))
p9#CDI


findings_0.05=apply((q_sum_CDI<0.05),2,as.numeric)
findings_0.05=as.data.frame(findings_0.05)
rownames(findings_0.05)=rownames(q_sum_CDI)

f1=function(i){
  t=0.05
  inter=rep(0,length(q_sum_CDI[,i]))
  inter[intersect(which(q_sum_CDI[,i]<t),which(q_sum_CDI[,i+k]<t))]=1
  inter
}

consist_0.05=foreach(i=1:k, .combine=cbind) %do% f1(i)
consist_0.05=as.data.frame(consist_0.05)
colnames(consist_0.05)=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
rownames(consist_0.05)=colnames(count)


###IGT###
###IGT_FULL###
count=IGT_KarlssonFH[[1]]
group=IGT_KarlssonFH[[2]]

adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC=NULL
}else{
  q_ANCOMBC <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM=NULL
}else{
  q_ANCOM=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda=NULL
}else{
  q_mbDecoda=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA=NULL
}else{
  q_LinDA=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM=NULL
}else{
  q_LOCOM=t(out_LOCOM[["q.otu"]])
}

###IGT_SUB###
count=IGT_KarlssonFH_sub[[1]]
group=IGT_KarlssonFH_sub[[2]]
adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC1=NULL
}else{
  q_ANCOMBC1 <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM1=NULL
}else{
  q_ANCOM1=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda1=NULL
}else{
  q_mbDecoda1=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA1=NULL
}else{
  q_LinDA1=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM1=NULL
}else{
  q_LOCOM1=t(out_LOCOM[["q.otu"]])
}

q_sum_IGT=data.frame(q_ANCOMBC,q_ANCOM,q_LinDA,q_LOCOM,q_mbDecoda,q_ANCOMBC1,q_ANCOM1,q_LinDA1,q_LOCOM1,q_mbDecoda1)
k=length(q_sum_IGT)/2


alpha=0.05
f=function(i){
  inter=length(intersect(which(q_sum_IGT[,i]<alpha),which(q_sum_IGT[,i+k]<alpha)))
  consist=inter/length(union(which(q_sum_IGT[,i]<alpha),which(q_sum_IGT[,i+k]<alpha)))
  consist
}


consist=foreach(i=1:k, .combine=rbind) %do% f(i)
consist[is.na(consist)]=0
method=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
out_IGT=data.frame(method,consist)


p10=ggplot(data=out_IGT,aes(method,consist,color=method,fill=method))+
  geom_bar(stat="identity",position=position_dodge(0.75),width = 0.5, color = "transparent")+ ylim(0,0.4)+
  theme_bw()+scale_fill_npg()+
  theme(plot.title = element_text(hjust=0.5),legend.position = "none", axis.text.x = element_blank())+
  labs(x="method",y="Jaccard",title = expression("IGT"))
p10#IGT



###ME/CFS###
###ME/CFS_FULL###
count=ME_NagySzakalD[[1]]
group=ME_NagySzakalD[[2]]

adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC=NULL
}else{
  q_ANCOMBC <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM=NULL
}else{
  q_ANCOM=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda=NULL
}else{
  q_mbDecoda=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA=NULL
}else{
  q_LinDA=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM=NULL
}else{
  q_LOCOM=t(out_LOCOM[["q.otu"]])
}

###ME/CFS_SUB###
count=ME_NagySzakalD_sub[[1]]
group=ME_NagySzakalD_sub[[2]]
adjust="BH"
Gamma=NULL
meta=data.frame(Sample.ID=rownames(count), group)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)

suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC1=NULL
}else{
  q_ANCOMBC1 <- out_ANCOMBC[["res"]][["q_val"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, zero_cut =1), silent = TRUE))
if (inherits(out_ANCOM, "try-error")) {
  q_ANCOM1=NULL
}else{
  q_ANCOM1=out_ANCOM[["results"]][["final"]][["log2FC.qval"]]
}

#mbDecoda
suppressWarnings(out_mbDecoda <- try(mbDecoda(count, x=group, Gamma=Gamma, W=NULL, prev.cut = 0, adjust=adjust)
                                  , silent = TRUE))
if (inherits(out_mbDecoda, "try-error")) {
  q_mbDecoda1=NULL
}else{
  q_mbDecoda1=out_mbDecoda[["Bias_correct"]][["DAA"]][["q.val"]]
}

#LinDA
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA1=NULL
}else{
  q_LinDA1=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM1=NULL
}else{
  q_LOCOM1=t(out_LOCOM[["q.otu"]])
}

q_sum_ME=data.frame(q_ANCOMBC,q_ANCOM,q_LinDA,q_LOCOM,q_mbDecoda,q_ANCOMBC1,q_ANCOM1,q_LinDA1,q_LOCOM1,q_mbDecoda1)
k=length(q_sum_ME)/2


alpha=0.05
f=function(i){
  inter=length(intersect(which(q_sum_ME[,i]<alpha),which(q_sum_ME[,i+k]<alpha)))
  consist=inter/length(union(which(q_sum_ME[,i]<alpha),which(q_sum_ME[,i+k]<alpha)))
  consist
}


consist=foreach(i=1:k, .combine=rbind) %do% f(i)
consist[is.na(consist)]=0
method=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
out_ME=data.frame(method,consist)


p11=ggplot(data=out_ME,aes(method,consist,color=method,fill=method))+
  geom_bar(stat="identity",position=position_dodge(0.75),width = 0.5, color = "transparent")+ ylim(0,0.4)+
  theme_bw()+scale_fill_npg()+
  theme(plot.title = element_text(hjust=0.5),legend.position = "none", axis.text.x = element_blank())+
  labs(x="method",y="Jaccard",title = expression("ME/CFS"))
p11#ME/CFS




#############numbers of finding
num=data.frame(method=rep(c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda"),2),
               findings=c(colSums(q_sum_CDI[,1:k]<alpha),colSums(q_sum_CDI[,-c(1:k)]<alpha)),
               dataset=rep(c("Schubert et al.","Vincent et al."),each=5))

num$findings=as.numeric(num$findings)
num$method=factor(num$method, levels =c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda"))
p12=ggplot(num, aes(y=findings, x=dataset,col=method)) +
  geom_point(size = 4)+
  geom_line(aes(group = method))+
  geom_vline(xintercept=c(1,2),linetype = "dashed")+
  scale_color_npg()+
  theme_classic()+labs(x=NULL,y= 'no. of detected taxa', col= "Method")
p12



p=plot_list(p1,p5,p2,p9,p3,p10,p4,p11,ncol =2,nrow = 4)|plot_list(p12,as.patchwork(p8),nrow =2)
p
