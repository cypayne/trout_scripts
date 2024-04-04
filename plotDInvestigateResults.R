# looking through Dinvestigate results:
# D has large variance when applied to small genomic windows and is a poor
# estimator of the amount of introgression
# --> see really high D values across genome, this tracks
# f_d: can be used to locate regions with introgression between P2 and P3
#      but can take arbitrarily large neg values when there's excess allele
#      sharing between P1 and P3
# f_dm: addresses this issue, equally quantifies shared variation between
#       P2 and P3 (positive) or P1 and P3 (negative)
# d_f:  uses correlation of allele frequencies (like f_d,f_dm) and also uses
#       genetic distance --> gives distance fraction

# Arlee genome info: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002163495.1/
# mitogenome: NC_001717.1
# cyp 12/2023

# read in the results with 50 SNP windows and a step of 25 SNPs
bigStep_all_clark_agua <- read.table("Dsuite_outs/Dinvestigate_windows/henshawii_1_clarkii_aguabonita_localFstats_clarkii-mykiss_trios_50_25.txt",as.is=T,header=T)
bigStep_all_clark_gaird <- read.table("Dsuite_outs/Dinvestigate_windows/henshawii_1_clarkii_gairdnerii_localFstats_clarkii-mykiss_trios_50_25.txt",as.is=T,header=T)
bigStep_all_clark_irid1 <- read.table("Dsuite_outs/Dinvestigate_windows/henshawii_1_clarkii_irideus_1_localFstats_clarkii-mykiss_trios_50_25.txt",as.is=T,header=T)
bigStep_all_clark_irid3 <- read.table("Dsuite_outs/Dinvestigate_windows/henshawii_1_clarkii_irideus_3_localFstats_clarkii-mykiss_trios_50_25.txt",as.is=T,header=T)

bigStep_all <- bigStep_all_clark_irid1
# match chr names to nums
chr_nums_table <- read.table("Arlee-chr_names2nums.txt",as.is=T,header=T)
bigStep_all_nums <- merge(bigStep_all,chr_nums_table,by.x='chr',by.y='chr_name')
bigStep_main <- subset(bigStep_all_nums, chr_num %in% (1:32))
# grab first 32 chr names
#chrom_names <- unique(bigStep_all$chr)[1:32]

par(mfrow = c(6, 5))
for( c in (1:32) ) {
  # subset chromosome
  bigStep <- subset(bigStep_main, chr_num==c)
  
  if(nrow(bigStep)==0) { next }
  
  # plot D in windows:
  #plot(bigStep$windowStart/1000000, bigStep$D,type="l",xlab="coordinate (Mb)",ylab="D (ABBA-BABA)")

  # plot f_dM in windows along scaffold 18:
  plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM")
}

## plot specific regions of interest
bigStep <- subset(bigStep_all, chr=="NC_048579.1")
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlim=c(23,23.25),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.8))

bigStep <- subset(bigStep_all, chr=="NC_048573.1")
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlim=c(13.1,13.6),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.8))

bigStep <- subset(bigStep_all, chr=="NC_048592.1")
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlim=c(36.2,36.5),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.8))

bigStep <- subset(bigStep_all, chr=="NC_048593.1")
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlim=c(38.1,38.4),xlab="scaffold 18 coordinate (Mb)",ylab="f_dM",ylim=c(-0.2,0.8))

  
## subset regions with high positive (i.e. sharing between P2 and P3) d_fM and d_f for each of the trios
subset(bigStep_all_clark_agua, f_dM > 0.7)
subset(bigStep_all_clark_gaird, f_dM > 0.7)
subset(bigStep_all_clark_irid1, f_dM > 0.7)
subset(bigStep_all_clark_irid3, f_dM > 0.7)


## collect the chromosomes with replicate f_dM > 0.7 between clarkii-irideus1 and clarkii-irideus3
clark_irid1_chrs <- unique(subset(bigStep_all_clark_irid1, f_dM > 0.65)$chr)
clark_irid3_chrs <- unique(subset(bigStep_all_clark_irid3, f_dM > 0.65)$chr)
irid1_irid3_chr_overlap <- intersect(clark_irid1_chrs,clark_irid3_chrs)

# look at the f_dM profiles across trios for both for these chromosomes
for( c in irid1_irid3_chr_overlap) {
  par(mfrow = c(4,1))
  # subset chromosome
  bigStep <- subset(bigStep_all_clark_irid1, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_1")
  
  bigStep <- subset(bigStep_all_clark_irid3, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_3")
  
  bigStep <- subset(bigStep_all_clark_gaird, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - gairdnerii")
  
  bigStep <- subset(bigStep_all_clark_agua, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - aguabonita")
  
}

# look at the f_dM profiles across trios for a specific region on a chromosome
c <- "NC_048579.1" # chr 15
start_Mbp <- 23
end_Mbp   <- 23.5

subset(bigStep_all_clark_agua, chr==c & windowStart > start_Mbp*1000000 & windowEnd < end_Mbp*1000000)

par(mfrow = c(4,1))
# subset chromosome
bigStep <- subset(bigStep_all_clark_irid1, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_1")
  
bigStep <- subset(bigStep_all_clark_irid3, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_3")

bigStep <- subset(bigStep_all_clark_gaird, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - gairdnerii")

bigStep <- subset(bigStep_all_clark_agua, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - aguabonita")


# look at the d_f profiles across trios for these chromosomes
for( c in irid1_irid3_chr_overlap) {
  par(mfrow = c(4,1))
  # subset chromosome
  bigStep <- subset(bigStep_all_clark_irid1, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_d,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,1),ylab="f_dM",main="henshawii_1 - clarkii - irideus_1")
  
  bigStep <- subset(bigStep_all_clark_irid3, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_d,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,1),ylab="f_dM",main="henshawii_1 - clarkii - irideus_3")
  
  bigStep <- subset(bigStep_all_clark_gaird, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_d,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,1),ylab="f_dM",main="henshawii_1 - clarkii - gairdnerii")
  
  bigStep <- subset(bigStep_all_clark_agua, chr==c)
  plot(bigStep$windowStart/1000000, bigStep$f_d,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,1),ylab="f_dM",main="henshawii_1 - clarkii - aguabonita")
  
}

par(mfrow = c(4,1))
# subset chromosome
bigStep <- subset(bigStep_all_clark_irid1, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_1")

bigStep <- subset(bigStep_all_clark_irid3, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_3")

bigStep <- subset(bigStep_all_clark_gaird, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - gairdnerii")

bigStep <- subset(bigStep_all_clark_agua, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),xlim=c(start_Mbp,end_Mbp),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - aguabonita")


## plot all trios for specific chromosomes
c <- "NC_048569.1" # chr 5
c <- "NC_048579.1" # chr 15

par(mfrow = c(4,1))
# subset chromosome
bigStep <- subset(bigStep_all_clark_irid1, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_1")

bigStep <- subset(bigStep_all_clark_irid3, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - irideus_3")

bigStep <- subset(bigStep_all_clark_gaird, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - gairdnerii")

bigStep <- subset(bigStep_all_clark_agua, chr==c)
plot(bigStep$windowStart/1000000, bigStep$f_dM,type="l",xlab=paste0("chr ", c, " coordinate (Mb)"),ylim=c(-0.4,0.8),ylab="f_dM",main="henshawii_1 - clarkii - aguabonita")
