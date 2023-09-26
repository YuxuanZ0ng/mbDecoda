library(stringr)
library(phyloseq)

###########################CDI#############################
otutable=read.table(file = "cdi_schubert.otu_table.txt",sep="\t",header=T)
metadata=read.table(file = "cdi_schubert.metadata.txt",sep="\t",header=T)
otutable1=read.table(file = "cdi_vincent.otu_table.txt",sep="\t",header=T)
metadata1=read.table(file = "cdi_vincent.metadata.txt",sep="\t",header=T)

tax <- as.data.frame(strsplit(otutable[,1],";"))
tax1 <- as.data.frame(strsplit(otutable1[,1],";"))
tax <- t(tax)
tax1 <- t(tax1)
rownames(tax) <- otutable[,1]
rownames(tax1) <- otutable1[,1]
rownames(otutable) <- otutable[,1]
rownames(otutable1) <- otutable1[,1]


otutable <- otutable[,-1]
otutable1 <- otutable1[,-1]


case.id <- metadata[which(metadata$DiseaseState=="CDI"),1]
control.id <- metadata[which(metadata$DiseaseState=="H"),1]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]


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






case.id1 <- metadata1[which(metadata1$DiseaseState=="CDI"),1]
control.id1 <- metadata1[which(metadata1$DiseaseState=="H"),1]
data.id1 <- c(case.id1,control.id1)
group1 = c(rep(1,length(case.id1)),rep(0,length(control.id1)))
group1 <- group1[which(data.id1%in%colnames(otutable1))]
# substr(colnames(otutable),2,8)
otutable1 <- otutable1[,match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))]]

tax1 <- tax1[,-(7:8)]
taxname1 <- phyloseq(tax_table(tax1))
taxa_names(taxname1) <- rownames(otutable1)

otutable1=t(as(otutable1,"matrix"))



phy_filter <- phyloseq(tax_table(taxname1), otu_table((otutable1),taxa_are_rows=F))
phy.glom <- tax_glom(phy_filter, taxrank="ta6")#, NArm=F
otu.data.f <- phy.glom@otu_table
datasize <- nrow(otu.data.f)

prevalence = apply(as(otu.data.f, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = prevalence>  0.05

otu.data.f1=prune_taxa(keepOTUs, otu.data.f)
count1=otu.data.f1@.Data



tax_g <- as.data.frame(strsplit(colnames(count),";"))
tax1_g <- as.data.frame(strsplit(colnames(count1),";"))
tax_g <- t(tax_g)[,-(7:8)]
tax1_g <- t(tax1_g)[,-(7:8)]
colnames(count)=apply(tax_g, 1, paste, collapse = ";")
colnames(count1)=apply(tax1_g, 1, paste, collapse = ";")


common <- intersect(colnames(count), colnames(count1))
count <- count[,match(common,colnames(count))]
count1 <- count1[,match(common,colnames(count1))]


CDI_schubert=list(count,group)
CDI_vincent=list(count1,group1)
saveRDS(CDI_schubert,file="CDI_schubert.RDS")
saveRDS(CDI_vincent,file="CDI_vincent.RDS")







###########################CRC#############################
otutable=read.table(file = "crc_baxter.otu_table.txt",sep="\t",header=T)
metadata=read.table(file = "crc_baxter.metadata.txt",sep="\t",header=T)
otutable1=read.table(file = "crc_zackular.otu_table.txt",sep="\t",header=T)
metadata1=read.table(file = "crc_zackular.metadata.txt",sep="\t",header=T, fileEncoding="ISO-8859-1")



tax <- as.data.frame(strsplit(otutable[,1],";"))
tax1 <- as.data.frame(strsplit(otutable1[,1],";"))
tax <- t(tax)
tax1 <- t(tax1)
rownames(tax) <- otutable[,1]
rownames(tax1) <- otutable1[,1]
rownames(otutable) <- otutable[,1]
rownames(otutable1) <- otutable1[,1]


otutable <- otutable[,-1]
colnames(otutable) <- as.data.frame(strsplit(colnames(otutable),"X"))[2,]
otutable1 <- otutable1[,-1]


case.id <- metadata[which(metadata$DiseaseState=="CRC"),1]
control.id <- metadata[which(metadata$DiseaseState=="H"),1]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]


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



case.id1 <- metadata1[which(metadata1$DiseaseState=="CRC"),1]
control.id1 <- metadata1[which(metadata1$DiseaseState=="H"),1]
data.id1 <- c(case.id1,control.id1)
data.id1 <- gsub("-", ".", data.id1)
group1 = c(rep(1,length(case.id1)),rep(0,length(control.id1)))
group1 <- group1[which(data.id1%in%colnames(otutable1))]
# substr(colnames(otutable),2,8)
otutable1 <- otutable1[,match(data.id1,colnames(otutable1))[which(data.id1%in%colnames(otutable1))]]

tax1 <- tax1[,-(7:8)]
taxname1 <- phyloseq(tax_table(tax1))
taxa_names(taxname1) <- rownames(otutable1)

otutable1=t(as(otutable1,"matrix"))



phy_filter <- phyloseq(tax_table(taxname1), otu_table((otutable1),taxa_are_rows=F))
phy.glom <- tax_glom(phy_filter, taxrank="ta6")#, NArm=F
otu.data.f <- phy.glom@otu_table
datasize <- nrow(otu.data.f)

prevalence = apply(as(otu.data.f, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = prevalence>  0.05

otu.data.f1=prune_taxa(keepOTUs, otu.data.f)
count1=otu.data.f1@.Data



tax_g <- as.data.frame(strsplit(colnames(count),";"))
tax1_g <- as.data.frame(strsplit(colnames(count1),";"))
tax_g <- t(tax_g)[,-(7:8)]
tax1_g <- t(tax1_g)[,-(7:8)]
colnames(count)=apply(tax_g, 1, paste, collapse = ";")
colnames(count1)=apply(tax1_g, 1, paste, collapse = ";")


common <- intersect(colnames(count), colnames(count1))
count <- count[,match(common,colnames(count))]
count1 <- count1[,match(common,colnames(count1))]



CRC_baxter=list(count,group)
CRC_zackular=list(count1,group1)
saveRDS(CRC_baxter,file="CRC_baxter.RDS")
saveRDS(CRC_zackular,file="CRC_zackular.RDS")






###########################IGT###################
library(curatedMetagenomicData)
library(dplyr)
library(mia)
KarlssonFH_2013=curatedMetagenomicData("KarlssonFH_2013.relative_abundance", dryrun = FALSE, counts = TRUE)%>%
  mergeData()
altExps(KarlssonFH_2013) <-
  splitByRanks(KarlssonFH_2013)
KarlssonFH_2013_genus=altExp(KarlssonFH_2013,"genus")
otutable=KarlssonFH_2013_genus@assays@data@listData[["relative_abundance"]]
metadata=as.data.frame(KarlssonFH_2013_genus@colData@listData)
rownames(metadata)=KarlssonFH_2013_genus@colData@rownames
#metadata$antibiotics_current_use=="no"&
case.id <- rownames(metadata)[which(metadata$study_condition=="IGT")]
control.id <- rownames(metadata)[which(metadata$study_condition=="control")]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]
otutable=t(as(otutable,"matrix"))


set.seed(12345)
#index=rbinom(nrow(otutable), 1, 1/3)sample(1:nrow(otutable),100)JieZ_2017
index=sample(1:nrow(otutable),nrow(otutable)/2)
otutable1=otutable[index,]
group1=group[index]

datasize <- nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  0.05)

datasize1 <- nrow(otutable1)
prevalence1 = apply(as(otutable1, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize1)

keepOTUs1 = which(prevalence1>  0.05)
common=intersect(keepOTUs,keepOTUs1)

otutable = otutable[,common]
otutable1 = otutable1[,common]

IGT_KarlssonFH=list(otutable,group)
IGT_KarlssonFH_sub=list(otutable1,group1)
saveRDS(IGT_KarlssonFH,file="IGT_KarlssonFH.RDS")
saveRDS(IGT_KarlssonFH_sub,file="IGT_KarlssonFH_sub.RDS")





###########################ME/CFS###################
NagySzakalD_2017=curatedMetagenomicData("NagySzakalD_2017.relative_abundance", dryrun = FALSE, counts = TRUE)%>%
  mergeData()
altExps(NagySzakalD_2017) <-
  splitByRanks(NagySzakalD_2017)
NagySzakalD_2017_genus=altExp(NagySzakalD_2017,"genus")
otutable=NagySzakalD_2017_genus@assays@data@listData[["relative_abundance"]]
metadata=as.data.frame(NagySzakalD_2017_genus@colData@listData)
rownames(metadata)=NagySzakalD_2017_genus@colData@rownames
case.id <- rownames(metadata)[which(metadata$study_condition=="ME/CFS")]
control.id <- rownames(metadata)[which(metadata$study_condition=="control")]
data.id <- c(case.id,control.id)
group = c(rep(1,length(case.id)),rep(0,length(control.id)))
group <- group[which(data.id%in%colnames(otutable))]
otutable <- otutable[,match(data.id,colnames(otutable))[which(data.id%in%colnames(otutable))]]
otutable=t(as(otutable,"matrix"))


set.seed(12345)
index=sample(1:nrow(otutable),nrow(otutable)/2)
otutable1=otutable[index,]
group1=group[index]

datasize <- nrow(otutable)
prevalence = apply(as(otutable, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize)

keepOTUs = which(prevalence>  0.05)

datasize1 <- nrow(otutable1)
prevalence1 = apply(as(otutable1, "matrix"), 2, function(x) {
  return(sum(x > 0))
})/(datasize1)

keepOTUs1 = which(prevalence1>  0.05)
common=intersect(keepOTUs,keepOTUs1)

otutable = otutable[,common]
otutable1 = otutable1[,common]

ME_NagySzakalD=list(otutable,group)
ME_NagySzakalD_sub=list(otutable1,group1)
saveRDS(ME_NagySzakalD,file="ME_NagySzakalD.RDS")
saveRDS(ME_NagySzakalD_sub,file="ME_NagySzakalD_sub.RDS")

