## NGSadmix outputs --> plots of selected K structure
## working with NGSadmix outputs from mega-post-bcf-exploratory-snakeflows pipeline
## output files for each K and rep are in folders called: K_#_rep_#
## R code adapted from Bay lab tutorial: 
## https://baylab.github.io/MarineGenomics/week-9-population-structure-using-ngsadmix.html
## cyp II-2024

ngsadmix_dir <- "ngsadmix/filt_snps05-thin100_0-maf05" # set NGSadmix outputs directory
N_K <- 10    # set number of K run
N_reps <- 4  # set number of reps run

## choose K

# pull all log files
log_files <- list.files(ngsadmix_dir, pattern = ".log", full.names = T, recursive=T)

# read in all logs
all_logs <- lapply(1:length(log_files), FUN = function(i) readLines(log_files[i]))

# make list of the line that starts with "best like=" from all logs, just target 'b'
library(stringr)
bestlikes_str_list <- sapply(1:length(log_files), FUN= function(x) all_logs[[x]][which(str_sub(all_logs[[x]], 1, 1) == 'b')])

# make dataframe with 1:N_K and N_reps to add likelihood values
loglikes <- data.frame(K = rep(1:N_K, each=N_reps))

# add the log likelihood (first number in the string)
loglikes$loglike<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", bestlikes_str_list) ))

# calculate delta K and probability
# choose the K with the highest value
# ex.
#   1         2         3         4         5 
# Inf       Inf  41.50002  29.32132 394.18113 
# in this case, it's 5 - will need to run more K so that we make sure best value is not >5
# ex. run 2
#   1           2           3           4           5           6           7           8           9          10 
# Inf   282.87382         Inf    41.50002    29.32132   394.18113         Inf 15828.78428   282.91030   202.86797 
# K=8 should be chosen
tapply(loglikes$loglike, loglikes$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

# chosen K
chosen_K   <- 8
chosen_rep <- 1

## plot

# read in the qopt ancestry proportions 
q <- read.table(file.path(ngsadmix_dir,paste0("K_",chosen_K,"_rep_",chosen_rep),"output.qopt"))

# read the sample_list.txt file; the qopt file is ordered by the sample order in this file
sample_list <- read.table(file.path(ngsadmix_dir,"sample_list.txt"))
colnames(sample_list) <- c("sample_ID")

# read metadata
metadata <- read.csv("refs/CCGP-trout_sample-coordinate-match.csv",header=T)

# merge sample list with metadata (in order of sample list)
sample_info <- merge(sample_list,metadata,by="sample_ID")

#order by population
ord <- order(sample_info$geo_rank_by_species)

# produce barplot
pdf("Figures/ccgp-trout_admix_plot_map-matched-colors_K7.pdf",width=20,height=8)
  par(mfrow=c(2,1))
  barplot(t(q)[,ord],
        col=hue_pal()(8),
        names=sample_info$new_sample_name[ord],
        las=2,
        space=0,
        border=NA,
        cex.names=0.8,
        xlab="",
        ylab=paste0("admixture proportions for K=",chosen_K))
dev.off()

par(mfrow=c(2,1))
barplot(t(q)[,ord],
        col=1:chosen_K,
        names=sample_info$plotting_locality[ord],
        las=2,
        space=0,
        border=NA,
        cex.names=0.8,
        xlab="",
        ylab=paste0("admixture proportions for K=",chosen_K))

barplot(t(q)[,ord],
        col=1:chosen_K,
        names=sample_info$subspecies_name[ord],
        las=2,
        space=0,
        border=NA,
        cex.names=0.8,
        xlab="",
        ylab=paste0("admixture proportions for K=",chosen_K))


# Snippet from Nina's code
admix.id = as.data.frame(cbind(pop, admix))
names(admix.id) = c("pop","q1","q2")

pdf(paste0(basedir, "/results/NGSadmix_LDpruned_K2_plot.pdf"))
plot = barplot(t(as.matrix(subset(admix.id, select=q1:q2))), col=c("firebrick","royalblue"), border=NA)
dev.off() 


## plot all reps 
total_reps <- 4

# prepare metadata
# read the sample_list.txt file; the qopt file is ordered by the sample order in this file
sample_list <- read.table(file.path(ngsadmix_dir,"sample_list.txt"))
colnames(sample_list) <- c("sample_ID")
# read metadata
metadata <- read.csv("refs/CCGP-trout_sample-coordinates_correctly-geo-ranked.csv",header=T)
# merge sample list with metadata (in order of sample list)
sample_info <- merge(sample_list,metadata,by="sample_ID")


# trout palette order: lahontan,  mccloud, steelhead, redband, coastal cutts, paiute, redband, coastal cutts, kern
trout_palette <- c("#D3928C", "#BBD3D1", "#7794A1","#546A7C", "#8E7068", "#253333", "#BFBDD1", "#6D7046")
for(k in c(2,3,4,5,6,7,8)) {
  ## plot
  pdf(paste0("Figures/all-ccgp_admix_plot_map-matched-colors_K",k,".pdf"),width=14,height=10)
  par(mfrow=c(total_reps+2,1), mar=c(1,4,1,1))
  for(rep in c(1:total_reps)) {
    # read in the qopt ancestry proportions 
    q <- read.table(file.path(ngsadmix_dir,paste0("K_",k,"_rep_",rep),"output.qopt"))
    q$sample_ID <- sample_list$sample_ID
    q$replicate <- rep
    #if(rep==1) { q_all_reps = q }
    #else { q_all_reps <- rbind(q_all_reps,q) 
    
    #order by geography
    ord <- order(sample_info$overall_rank)
    
    if(rep==total_reps) {
      barplot(t(q)[,ord],
              #        col=hue_pal()(chosen_K),
              col=trout_palette,
              names=sample_info$new_sample_name[ord],
              las=2,
              space=0,
              border=NA,
              cex.names=1,
              xlab="",
              ylab=paste0("admixture proportions for K=",k)) 
    }
    else {
      barplot(t(q)[,ord],
              #        col=hue_pal()(chosen_K),
              col=trout_palette,
              las=2,
              space=0,
              border=NA,
              xlab="",
              ylab=paste0("admixture proportions for K=",k)) 
    }
  }
  dev.off()
}





### Map with ancestry pie charts
## to add piecharts of ancestry to a map
## use scatterpie? https://cran.r-project.org/web/packages/scatterpie/vignettes/scatterpie.html#scatter-pie-plot
## data should be in this format
##          long        lat region          A        B        C        D
## 1  -56.047565  12.665926      1 2.13121969 8.663359 3.928711 8.676792
## 2  -23.017749  -1.427338      2 0.25688371 1.403569 1.375096 4.945092
## 4    7.050839  68.430114      3 0.24669188 0.524395 3.189978 5.138863

## format data
## create dataframe of sample names and their ancestry proportions

# read the sample_list.txt file; the qopt file is ordered by the sample order in this file
ngsadmix_dir <- "ngsadmix/filt_snps05-thin100_0-maf05" # set NGSadmix outputs directory
sample_list <- read.table(file.path(ngsadmix_dir,"sample_list.txt"))
colnames(sample_list) <- c("sample_ID")

# read pop and coordinate metadata
metadata <- read.csv("refs/CCGP-trout_sample-coordinate-match.csv",header=T)

# merge sample list with metadata (in order of sample list)
sample_info <- merge(sample_list,metadata,by="sample_ID")

# read in the qopt ancestry proportions 
chosen_K   <- 8
chosen_rep <- 1
admix_props <- read.table(file.path(ngsadmix_dir,paste0("K_",chosen_K,"_rep_",chosen_rep),"output.qopt"))
names(admix_props) = paste0("q",1:chosen_K)

# merge sample-coordinate metadata with ancestry proportions data
sample_coords_admix <- as.data.frame(cbind(sample_info, admix_props))

# remove the hatchery fish
sample_coords_admix_nonNA_nohatchery <- subset(sample_coords_admix_nonNA, subspecies_name != "hatchery")

library(maps)
library(ggplot2)
library(scatterpie)
#https://rdrr.io/cran/maps/man/map.html
#ca_map <- map('county', 'california')

counties <- map_data("county")
ca_county <- subset(counties, region == "california")
ca_base <- ggplot(data = ca_county, mapping = aes(x = long, y = lat, group = region)) + 
  geom_polygon(color = "white", fill = "lightgrey") + coord_fixed(1.3) +
  geom_path(data=ca_rivers.sf_Rivers_coords.df, aes(X, Y, group=L1, fill=NULL), 
            color="#2A788EFF", alpha=0.7) 
aesthetic_map <-
  ca_base + theme_bw() +
  geom_scatterpie(aes(x=longitude, y=latitude, group=geo_group), data=sample_coords_admix_nonNA_nohatchery, cols=paste0("q",1:chosen_K), color=NA) +
  geom_label(data=sample_coords_admix_nonNA_nohatchery, aes(x=longitude, y=latitude, group=geo_group, label=geo_group), 
             nudge_x=0, nudge_y=-0, label.size=0.2, size=2,  
             fontface = "bold.italic", label.r=unit(0.20, "lines")) 
#          scale_fill_manual(value=trout_palette)
aesthetic_map

# add a scale bar and north arrow
# https://stackoverflow.com/questions/61809382/how-can-i-put-a-scalebar-and-a-north-arrow-on-the-map-ggplot
library(ggspatial)

ca_base +
  theme_bw() +
  geom_scatterpie(aes(x=longitude, y=latitude, group=region), data=sample_coords_admix_nonNA_nohatchery, cols=paste0("q",1:chosen_K), color=NA) +
  geom_label(data=sample_coords_admix_nonNA_nohatchery, aes(x=longitude, y=latitude, group=region, label=sample_ID), 
             nudge_x=0.4, nudge_y=-0.1, label.size=0.05,size=2, 
             fontface = "bold.italic", label.r=unit(0.20, "lines")) +
  ggspatial::annotation_scale(location = "tr", bar_cols = c("grey60", "white")) +
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "true",
    pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"))

ggsave(aesthetic_map, filename='Figures/ccgp-trout_admixture_map_aesthetic.pdf', device = cairo_pdf, width=6, height=10, bg="transparent")


# temp 
ggplot() + geom_scatterpie(aes(x=longitude, y=latitude, group=region), data=sample_coords_admix_nonNA, cols=paste0("q",1:chosen_K)) + coord_equal()


shapeUFs <- readOGR('.', 'BRUFE250GC_SIR')
shapeHid <- readOGR('.', 'PrincipaisRiosDoBrasil') 

ggplot(ca_rivers, aes(long, lat, group = group)) +
  geom_polygon(fill = 'gray90', color = 'black') +
  geom_path(data = shapeHid, color = 'steelblue2') +
  coord_map() + theme_void()









