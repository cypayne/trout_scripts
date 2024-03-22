## clustering with dendrograms and pca
## using outputs from snpgds on hpc
## see below for the snpgds code run on hpc

library(SNPRelate)
library(dplyr)
library(ggplot2)
library(ggrepel)

# set trout palette
# {"Tea rose (red)":"F6C7C0","Cadet gray":"96ACC0","Reseda green":"6F6B4B","Lavender (web)":"DDD4DF","Ash gray":"BFD8D4","Outer space":"324141","Cinereous":"88716C"}
trout_palette <- c("#96ACC0", "#88716C", "#F6C7C0", "#324141", "#DDD4DF","#BFD8D4","#6F6B4B")

# open gds file
genofile<-snpgdsOpen("snpgds_dendrograms/all-ccgp-trout-with-Chinook-outgroup.2Arlee.32chr-nums.vcf.gds")

#####
## code for creating dendrogram
## for reference, here is the (computationally expensive) code run on an hpc
## adapted from https://benbowlab.github.io/Benin-NGS/

library(gdsfmt)
library(SNPRelate)
library(ggplot2)

# convert vcf to gds format (for memory efficiency)
vcf.in <- "all-ccgp-trout-with-Chinook-outgroup.2Arlee.32chr-nums.vcf"
snpgdsVCF2GDS(vcf.in,paste0(vcf.in,".gds"),method ="biallelic.only")

# make relatedness matrix
genofile<-snpgdsOpen("all-ccgp-trout-with-Chinook-outgroup.2Arlee.32chr-nums.vcf.gds")
set.seed(100)
ibs.hc<-snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))
saveRDS(ibs.hc,"ibs_clusters.rds")
rv <- snpgdsCutTree(ibs.hc)

saveRDS(rv,"relatedness_tree.rds")
# plot relatedness dendrogram
pdf("snpgds_dendrogram.pdf",width=14, height=10)
plot(rv$dendrogram,main="ccgp trout SNP dendrogram")
dev.off()

# make dissimilarity matrix
dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, autosome.only=FALSE,remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)
saveRDS(dissMatrix,"dissimilarity_matrix.rds")
snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)
cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL, col.outlier="red", col.list=NULL, pch.outlier=4, pch.list=NULL,label.H=FALSE, label.Z=TRUE, verbose=TRUE)
saveRDS(cutTree,"dissimilarity_tree.rds")

# plot dissimilarity dendrogram
pdf("snpgds_dissimilarity_dendrogram.pdf",width=10,height=10)
snpgdsDrawTree(cutTree, main = "Phylogenetic Tree",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"), y.label.kinship=T,leaflab="perpendicular")
dev.off()
