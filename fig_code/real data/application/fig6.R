library(curatedMetagenomicData)
library(dplyr)
library(mia)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(phyloseq)
library(ANCOMBC)
library(MicrobiomeStat)
library(glmmTMB)
library(fastANCOM)
library(LOCOM)
library(UpSetR)
library(aplot)
library(patchwork)
library(doSNOW)
library(ggVennDiagram)
source("mbDecoda.R")



################  OB ############################
otutable=read.table(file = "ob_zupancic.otu_table.txt",sep="\t",header=T)
metadata=read.table(file = "ob_zupancic.metadata.txt",sep="\t",header=T)

tax <- as.data.frame(strsplit(otutable[,1],";"))
tax <- t(tax)
rownames(tax) <- otutable[,1]
rownames(otutable) <- otutable[,1]
otutable <- otutable[,-1]

case.id <- metadata[which(metadata$DiseaseState=="OB"),1]
control.id <- metadata[which(metadata$DiseaseState=="H"),1]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]
Gamma = metadata[match(colnames(otutable),metadata$X),]$sex_s
Gamma <- as.factor(Gamma) 


tax <- tax[,-(7:8)]
taxname <- phyloseq(tax_table(tax))
taxa_names(taxname) <- rownames(otutable)
otutable=t(as(otutable,"matrix"))


phy_filter <- phyloseq(tax_table(taxname), otu_table((otutable),taxa_are_rows=F))
phy.glom <- tax_glom(phy_filter, taxrank="ta6")#, NArm=F
otu.data.f <- phy.glom@otu_table
datasize <- nrow(otu.data.f)

prevalence = apply(as(otu.data.f, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = prevalence>  0.05

otu.data.f1=prune_taxa(keepOTUs, otu.data.f)
count=otu.data.f1@.Data


tax_g <- as.data.frame(strsplit(colnames(count),";"))
tax_g <- t(tax_g)[,-(7:8)]
colnames(count)=apply(tax_g, 1, paste, collapse = ";")



OB_zupancic=list(count,group,Gamma)



count=OB_zupancic[[1]]
group=OB_zupancic[[2]]
Gamma=OB_zupancic[[3]]
adjust="BH"
meta=data.frame(Sample.ID=rownames(count), group,x1=Gamma)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)


suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group+x1",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC=NULL
}else{
  q_ANCOMBC <- out_ANCOMBC[["res"]][["q_val"]][["group"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, Z=Gamma, zero_cut =1), silent = TRUE))
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
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group+x1',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA=NULL
}else{
  q_LinDA=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), C=Gamma, fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM=NULL
}else{
  q_LOCOM=t(out_LOCOM[["q.otu"]])
}

q_sum_OB=data.frame(q_ANCOMBC,q_ANCOM,q_LinDA,q_LOCOM,q_mbDecoda)
findings=apply((q_sum_OB<0.05),2,as.numeric)
findings=as.data.frame(findings)
colnames(findings)=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
rownames(findings)=colnames(count)


id <- list(
  ANCOMBC=which(findings$ANCOMBC==1),
  fastANCOM=which(findings$fastANCOM ==1),
  LinDA=which(findings$LinDA==1),
  LOCOM=which(findings$LOCOM ==1),
  mbDecoda=which(findings$mbDecoda==1)
)

p1=ggVennDiagram(id,set_size = 3,
                 category.names = c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda"),label = "count", 
                 label_color = "black",
                 label_alpha = 0,
                 edge_lty = "dashed", 
                 edge_size = 1)+ scale_fill_gradient(low="#1d629510",high = "#2a83a2",guide="none")+ 
  scale_x_continuous(expand = expansion(mult = .2))+
  theme(plot.title = element_text(hjust=0.5))+
  labs(title = expression(underline("OB")))
p1




diff=which(rowSums(findings)==1&findings$mbDecoda==1)
df=q_sum_OB[diff,]
f=function(i){
  ta=unlist(strsplit(i, "[;]"))
  ta=ta[nchar(ta) > 3]
  c(gsub("p__", "", ta[2]),gsub("c__", "",gsub("g__", "", ta[length(ta)])))
}
df=data.frame(foreach(i=rownames(df), .combine=rbind) %do% f(i),out_mbDecoda[["Bias_correct"]][["DAA"]][["delta"]][diff],q_mbDecoda[diff],"A")
colnames(df)=c("Phylum","Taxa","logFC","q","x")
df$Phylum_col[which(df$Phylum == "Firmicutes")] <- "#e7a40e"
df$Phylum_col[which(df$Phylum != "Firmicutes")] <- "#78bee5"
df$q[df$q<0.01]="**"
df$q[df$q<0.05&df$q>=0.01]="*"
df$q_col[which(df$q == "*" & df$logFC > 0)] <- "Postive effect(P<0.05)"
df$q_col[which(df$q == "**" & df$logFC > 0)] <- "Postive effect(P<0.01)"
df$q_col[which(df$q == "*" & df$logFC < 0)] <- "Negtive effect(P<0.05)"
df$q_col[which(df$q == "**" & df$logFC < 0)] <- "Negtive effect(P<0.01)"
df$Phylum <- factor(df$Phylum, levels = rev(unique(df$Phylum)))
df$Taxa <- factor(df$Taxa, levels = rev(unique(df$Taxa)))
df <- df[order(df$Phylum), ]


p2=ggplot(df)+
  # 0轴竖线：
  geom_hline(yintercept = 0, linewidth = 0.3)+
  # 线条???
  geom_segment(aes(y = 0, x = Taxa, yend = logFC, xend = Taxa, color = q_col),size = 1) +
  # 散点???
  geom_point(aes(Taxa, logFC, color = q_col)) +
  # 显著性：
  geom_text(aes(x= as.numeric(Taxa) + 0.1, y = logFC + 0.2, label = q, color = q_col), show.legend = F)+
  # 颜色???
  scale_color_manual(name = "",values = c("Postive effect(P<0.05)" = "#d55e00",
                                          "Postive effect(P<0.01)" = "#ffbd88",
                                          "Negtive effect(P<0.05)" = "#0072b2",
                                          "Negtive effect(P<0.01)" = "#7acfff"))+ 
  # 背景色：
  annotate("rect",
           xmin = c(0.3,4.5),
           xmax = c(4.5,6.7),
           ymin = -3.5, ymax = 3, alpha = 0.2, fill = rev(unique(df$Phylum_col)))+
  # 调整x轴拓宽：
  scale_y_continuous(expand = c(0,0))+
  xlab("")+
  ylab("log fold change")+
  theme_bw()+
  guides(color=F)+ 
  theme(axis.text.y = element_text(color = rev(df$Phylum_col)))+
  coord_flip()+
  labs(title= expression(underline("OB")))

annotation_custom2 <- function (grob, 
                                xmin = -Inf, 
                                xmax = Inf, 
                                ymin = -Inf, 
                                ymax = Inf, 
                                data) 
{
  layer(data = data, 
        stat = StatIdentity, 
        position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

p2=p2+annotation_custom2(
  data = df %>% filter(Phylum=="Proteobacteria"),
  grob = grid::textGrob(label="Proteobacteria",
                        rot=90,
                        gp=grid::gpar(fontsize= 9,col='#0072b2')),
  xmin=5.6,xmax=5.6,
  ymin=-3.35,ymax=-3.35)+
  annotation_custom2(
    data = df %>% filter(Phylum=="Firmicutes"),
    grob = grid::textGrob(label="Firmicutes",
                          rot=90,
                          gp=grid::gpar(fontsize= 9,col='#d55e00')),
    xmin=2.5,xmax=2.5,
    ymin=-3.35,ymax=-3.35)+
  theme(plot.title = element_text(hjust = 0.5))
p2


################  T2D ############################
metadata=read.csv(file = "T2D_meta.csv",header=T)
QinJ_2012=curatedMetagenomicData("QinJ_2012.relative_abundance", dryrun = FALSE, counts = TRUE)%>%
  mergeData()
altExps(QinJ_2012) <-
  splitByRanks(QinJ_2012)
QinJ_2012_genus=altExp(QinJ_2012,"genus")
otutable=QinJ_2012_genus@assays@data@listData[["relative_abundance"]]

case.id <- metadata[which(metadata$Diabetic..Y.or.N.=="Y"),1]
control.id <- metadata[which(metadata$Diabetic..Y.or.N.=="N"),1]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]
Gamma = metadata[match(colnames(otutable),metadata$Sample.ID),18]
Gamma <- as.factor(Gamma) 
otutable=t(as(otutable,"matrix"))


datasize <- nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  0.05)
otutable = otutable[,keepOTUs]
#Gamma=NULL

T2D_QinJ=list(otutable,group,Gamma)



count=T2D_QinJ[[1]]
group=T2D_QinJ[[2]]
Gamma=T2D_QinJ[[3]]
adjust="BH"
meta=data.frame(Sample.ID=rownames(count), group,x1=Gamma)
rownames(meta)=rownames(count)
OTU = otu_table(t(count), taxa_are_rows = TRUE)
META = sample_data(meta)
PHYSEQ = phyloseq(OTU, META)


suppressWarnings(out_ANCOMBC <- try(ancombc(phyloseq = PHYSEQ,formula ="group+x1",
                                            p_adj_method = adjust, zero_cut = 1.1, lib_cut = 1), 
                                    silent = TRUE))
if (inherits(out_ANCOMBC, "try-error")) {
  q_ANCOMBC=NULL
}else{
  q_ANCOMBC <- out_ANCOMBC[["res"]][["q_val"]][["group"]]
}

suppressWarnings(out_ANCOM <- try(fastANCOM(Y=as.matrix(count), x=group, Z=Gamma, zero_cut =1), silent = TRUE))
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
suppressWarnings(out_LinDA <- try(linda(t(count), meta,formula = '~group+x1',alpha = 0.05,p.adj.method = adjust)
                                  , silent = TRUE))
if (inherits(out_LinDA, "try-error")) {
  q_LinDA=NULL
}else{
  q_LinDA=out_LinDA[["output"]][["group"]][["padj"]]
}

# LOCOM
suppressWarnings(out_LOCOM <- try(locom(otu.table =as.matrix(count), Y = factor(group), C=Gamma, fdr.nominal = 0.05, prev.cut = 0, seed = 1, n.cores = 1)
                                  , silent = TRUE))
if (inherits(out_LOCOM, "try-error")) {
  q_LOCOM=NULL
}else{
  q_LOCOM=t(out_LOCOM[["q.otu"]])
}




q_sum_T2D=data.frame(q_ANCOMBC,q_ANCOM,q_LinDA,q_LOCOM,q_mbDecoda)
findings=apply((q_sum_T2D<0.05),2,as.numeric)
findings=as.data.frame(findings)
colnames(findings)=c("ANCOMBC","fastANCOM","LinDA","LOCOM","mbDecoda")
rownames(findings)=colnames(count)

queries = list(list(query = intersects,
                    params = list("mbDecoda"),
                    active = T,
                    color="#D55E00",
                    query.name = "ABC"))


p3=UpSetR::upset(findings,
      point.size = 1.8,
      line.size = 1,
      main.bar.color = "#2a83a2",
      sets.bar.color = c("#3b7960","#3b7960","#3b7960","#D55E00","#3b7960"),
      queries = queries)
p3


diff=which(rowSums(findings)==1&findings$mbDecoda==1)
df1=q_sum_T2D[diff,]
f=function(i){
  ta=unlist(strsplit(i, "[_]"))
  ta=ta[nchar(ta) > 3]
  c(ta[1],ta[length(ta)])
}
df1=data.frame(foreach(i=rownames(df1), .combine=rbind) %do% f(i),out_mbDecoda[["Bias_correct"]][["DAA"]][["delta"]][diff],q_mbDecoda[diff],"A",NA)
colnames(df1)=c("Phylum","Taxa","logFC","q","x","Phylum_col")
df1$Phylum_col[which(df1$Phylum == "Firmicutes")] <- "#e7a40e"
df1$Phylum_col[which(df1$Phylum != "Firmicutes")] <- "#2a83a2"
df1$q[df1$q<0.01]="**"
df1$q[df1$q<0.05&df1$q>=0.01]="*"
df1$q_col[which(df1$q == "*" & df1$logFC > 0)] <- "Postive effect(P<0.05)"
df1$q_col[which(df1$q == "**" & df1$logFC > 0)] <- "Postive effect(P<0.01)"
df1$q_col[which(df1$q == "*" & df1$logFC < 0)] <- "Negtive effect(P<0.05)"
df1$q_col[which(df1$q == "**" & df1$logFC < 0)] <- "Negtive effect(P<0.01)"
df1$Phylum <- factor(df1$Phylum, levels = rev(unique(df1$Phylum)))
df1$Taxa <- factor(df1$Taxa, levels = rev(unique(df1$Taxa)))
df1 <- df1[order(df1$Phylum), ]

p4=ggplot(df1)+
  # 0轴竖线：
  geom_hline(yintercept = 0, linewidth = 0.3)+
  # 线条???
  geom_segment(aes(y = 0, x = Taxa, yend = logFC, xend = Taxa, color = q_col),size = 1) +
  # 散点???
  geom_point(aes(Taxa, logFC, color = q_col)) +
  # 显著性：
  geom_text(aes(x= as.numeric(Taxa) + 0.1, y = logFC + 0.2, label = q, color = q_col), show.legend = F)+
  # 颜色???
  scale_color_manual(name = "",values = c("Postive effect(P<0.05)" = "#d55e00",
                                          "Postive effect(P<0.01)" = "#ffbd88",
                                          "Negtive effect(P<0.05)" = "#0072b2",
                                          "Negtive effect(P<0.01)" = "#7acfff"))+ 
  # 背景色：
  annotate("rect",
           xmin = c(0.3,4.5),
           xmax = c(4.5,5.8),
           ymin = -3.5, ymax = 3, alpha = 0.2, fill = rev(unique(df1$Phylum_col)))+
  # 调整x轴拓宽：
  scale_y_continuous(expand = c(0,0))+
  xlab("")+
  ylab("log fold change")+
  theme_bw()+
  guides(color=F)+
  labs(title= expression(underline("T2D")))+ 
  theme(axis.text.y = element_text(color = rev(df1$Phylum_col)))+
  coord_flip()


p4=p4+annotation_custom2(
  data = df1 %>% filter(Phylum=="Bacteroidota"),
  grob = grid::textGrob(label="Bacteroidota",
                        rot=90,
                        gp=grid::gpar(fontsize= 9,col='#2a83a2')),
  xmin=5.15,xmax=5.15,
  ymin=-3.35,ymax=-3.35)+
  annotation_custom2(
    data = df1 %>% filter(Phylum=="Firmicutes"),
    grob = grid::textGrob(label="Firmicutes",
                          rot=90,
                          gp=grid::gpar(fontsize= 9,col='#d55e00')),
    xmin=2.5,xmax=2.5,
    ymin=-3.35,ymax=-3.35)+
  theme(plot.title = element_text(hjust = 0.5))
p4

