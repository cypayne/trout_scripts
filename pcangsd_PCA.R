## PCA from pcangsd covariance matrix
## pulled from cluster

library(ggplot2)
library(ggrepel)


# set palette
trout_palette <- c("#546A7C","#BFBDD1","#8E7068","#253333","#BBD3D1","#F4AD8B","#D3928C","#6D7046","#7794A1")

cov_file_path <- "mykiss-only/mykiss-only_ccgp-plus-10-new-samps.cov"
sample_list_path <- "mykiss-only/sample_list.txt"
metadata_path <- "refs/CCGP-trout_sample-coordinate-match_plus-10-new-samples.csv"
# cov_file_path <- "ccgp-problem-trout/problem-ccgp-trout-pops_PCA.cov"
# metadata_path <- "metadatas/ccgp-problem-pop_pca-metadata.csv"

# read covariance matrix and metadata
cov <- as.matrix(read.table(cov_file_path))

# compute eigenvalues from covariance matrix
e<-eigen(cov)

# check variance explained by each component
pc_var <- e$values/sum(e$values)
pc_var

# simple plot
plot(e$vectors[,1:2])

par(mfrow=c(5,5))
for(i in c(1:5)) {
  for(j in c(1:5)) {
    if(i==j) next
    plot(e$vectors[,i:j])
  }
}

# make eigenvector dataframe
evectors <- as.data.frame(e$vectors)

# make a new dataframe with eigenvectors and metadata
sample_list <- read.table(sample_list_path)
# eigenvector rows are in the order of the list of samples
rownames(evectors)<-sample_list$V1

# merge eigenvector and metadata dataframes
metadata <- read.csv(metadata_path, header=T)
pca_data <- merge(metadata,evectors,by.x='sample_ID',by.y='row.names')


# grab the eigenvalues for labeling PC1 and PC2
x_pc <- 1
y_pc <- 2

X_PC_label <- round((pc_var[x_pc]*100),digits=2)
Y_PC_label <- round((pc_var[y_pc]*100),digits=2)
pca <- ggplot(data = pca_data, aes(x = V1, y = V2, color=subspecies_name, label=sample_ID)) +
  geom_point(aes(color = subspecies_name), alpha = 0.8, size=2) + 
  scale_color_manual(values=trout_palette) +
  geom_text_repel(aes(label = sample_ID), size = 3.5, max.overlaps=20) +
  theme_bw() +
  xlab(paste0("PC ", x_pc, " - ",X_PC_label,"%")) + ylab(paste0("PC ", y_pc, " - ",Y_PC_label,"%"))
pca

ggsave(pca, filename='Figures/ccgp-plus-10-new-samps_pc1-pc2.pdf', device = cairo_pdf, width=14, height=6, bg="transparent")
