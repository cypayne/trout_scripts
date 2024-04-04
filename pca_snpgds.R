## clustering with pca
## using outputs from snpgds on hpc
## see below for the snpgds code run on hpc

library(SNPRelate)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)


## TODO
# in the zoomed PCA replace x and y labels without the variation 
# PCA w just mykiss
# redo tree just with mykiss and no hatchery fish, see how trees compare
# hatchery locations are temporary stations, best to use strain name
#   a lot of them were sampled in ventura county
# sac river may not be as informative, populations are distinct along this 
#   huge river basin - use feather river
# add unplaced individuals

# set trout palette
# {"Tea rose (red)":"F6C7C0","Cadet gray":"96ACC0","Reseda green":"6F6B4B","Lavender (web)":"DDD4DF","Ash gray":"BFD8D4","Outer space":"324141","Cinereous":"88716C"}
trout_palette <- c("#96ACC0", "#88716C", "#F6C7C0", "#324141", "#DDD4DF","#BFD8D4","#6F6B4B","#8DA750","salmon","navy")

# open gds file
genofile<-snpgdsOpen("snpgds_dendrograms/all-ccgp-trout-with-Chinook-outgroup.2Arlee.32chr-nums.vcf.gds")

# create PCA
pca <- snpgdsPCA(genofile,remove.monosnp=T,maf=0.1,missing.rate=0.1) # #SNPs used: 12,784,578
plot(pca) # look at PC1 v PC2
plot(pca, 1:4) # look at plots for all combos between PC1 through PC4

# to overlay population information as colors
# load metadata
metadata <- read.table("refs/SETS-w-Chinook.plus.tsv",header=T)
metadata
# set subspecies of "unplaced#" samples to "unplaced" so they're binned
metadata$subspecies_name <- str_replace(metadata$subspecies_name, "unplaced\\d+", "unplaced")

# match the samples with metadata
pca_meta <- data.frame(sample.id = pca$sample.id,
                       subsp = factor(metadata$subspecies_name)[match(pca$sample.id, metadata$sample_ID)],
                       pop = factor(metadata$watershed_ID)[match(pca$sample.id, metadata$sample_ID)],
                       locality = factor(metadata$plotting_locality)[match(pca$sample.id, metadata$sample_ID)],
                       EV1 = pca$eigenvect[,1],    # the first eigenvector
                       EV2 = pca$eigenvect[,2],    # the second eigenvector
                       stringsAsFactors = FALSE)
pca_meta

PC1_label <- round((pca$varprop[1]*100),digits=2)
PC2_label <- round((pca$varprop[2]*100),digits=2)
all_pca_plot <-
  ggplot(data = pca_meta, aes(x = EV1, y = EV2, color=subsp, label=pop)) +
  geom_point(aes(color = subsp), alpha = 0.8, size=2) + 
  geom_text(hjust=-0.1, vjust=0) +
  scale_color_manual(values=trout_palette,name="subspecies") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12)) +
  xlab(paste0("PC 1 - ",PC1_label,"%")) + ylab(paste0("PC 2 - ",PC2_label,"%"))
all_pca_plot
ggsave(all_pca_plot,filename="Figures/PCA_ccgp_all.pdf",height=15,width=15,bg="transparent")


## rerun pca without chinook outgroup

just_trout_samples <- subset(metadata,subspecies_name!="chinook")$sample_ID
pca_just_trout <- snpgdsPCA(genofile,remove.monosnp=T,maf=0.1,missing.rate=0.1,sample.id=just_trout_samples) # # of SNPs: 12,804,807
# match the samples with metadata
pca_just_trout_meta <- data.frame(sample.id = pca_just_trout$sample.id,
                       subsp = factor(metadata$subspecies_name)[match(pca_just_trout$sample.id, metadata$sample_ID)],
                       pop = factor(metadata$watershed_ID)[match(pca_just_trout$sample.id, metadata$sample_ID)],
                       locality = factor(metadata$plotting_locality)[match(pca_just_trout$sample.id, metadata$sample_ID)],
                       EV1 = pca_just_trout$eigenvect[,1],    # the first eigenvector
                       EV2 = pca_just_trout$eigenvect[,2],    # the second eigenvector
                       stringsAsFactors = FALSE)
# grab the eigenvalues for labeling PC1 and PC2
PC1_label <- round((pca_just_trout$varprop[1]*100),digits=2)
PC2_label <- round((pca_just_trout$varprop[2]*100),digits=2)
ggplot(data = pca_just_trout_meta, aes(x = EV1, y = EV2, color=subsp, label=pop)) +
  geom_point(aes(color = subsp), alpha = 0.8, size=2) + 
  scale_color_manual(values=trout_palette) +
  #geom_text_repel(aes(label = pop), size = 3.5, max.overlaps=20) +
  theme_bw() +
  xlab(paste0("PC 1 - ",PC1_label,"%")) + ylab(paste0("PC 2 - ",PC2_label,"%"))


## rerun pca just on O mykiss
trout_palette_sub <- c("#96ACC0", "#6F6B4B", "#DDD4DF", "#324141","salmon")
just_mykiss_samples <- subset(metadata,species_name=="mykiss")$sample_ID
pca_just_mykiss <- snpgdsPCA(genofile,remove.monosnp=T,maf=0.1,missing.rate=0.1,sample.id=just_mykiss_samples) # # of SNPs: 5,326,135
# match the samples with metadata
pca_just_mykiss_meta <- data.frame(sample.id = pca_just_mykiss$sample.id,
                                  subsp = factor(metadata$subspecies_name)[match(pca_just_mykiss$sample.id, metadata$sample_ID)],
                                  pop = factor(metadata$watershed_ID)[match(pca_just_mykiss$sample.id, metadata$sample_ID)],
                                  locality = factor(metadata$plotting_locality)[match(pca_just_mykiss$sample.id, metadata$sample_ID)],
                                  long_name = factor(metadata$new_sample_name)[match(pca_just_mykiss$sample.id, metadata$sample_ID)],
                                  EV1 = pca_just_mykiss$eigenvect[,1],    # the first eigenvector
                                  EV2 = pca_just_mykiss$eigenvect[,2],    # the second eigenvector
                                  stringsAsFactors = FALSE)
# grab the eigenvalues for labeling PC1 and PC2
PC1_label <- round((pca_just_mykiss$varprop[1]*100),digits=2)
PC2_label <- round((pca_just_mykiss$varprop[2]*100),digits=2)
mykiss_pca_plot <- 
  ggplot(data = pca_just_mykiss_meta, aes(x = EV1, y = EV2, color=subsp, label=pop)) +
  geom_point(aes(color = subsp), alpha = 0.8, size=2) + 
  geom_text_repel(aes(label = locality), size = 3.5, max.overlaps=40) +
  scale_color_manual(values=trout_palette_sub,name="subspecies") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14)) +
  xlab(paste0("PC 1 - ",PC1_label,"%")) + ylab(paste0("PC 2 - ",PC2_label,"%"))
mykiss_pca_plot
ggsave(mykiss_pca_plot,filename="Figures/PCA_ccgp_Omykiss.pdf",height=10,width=18,bg="transparent")

plot(pca_just_mykiss, 1:4) # look at plots for all combos between PC1 through PC4


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
