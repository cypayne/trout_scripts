## Dsuite output exploration
## following analysis tutorial 
## https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data
library(dplyr)

# load Dmin results
infile <- "Dsuite_outs/broad-group_all-ccgp-nohatch-nobad-trout_dsuite/all-ccgp-nohatch-nobad-trout-with-Chinook-outgroup_broad-group_dsuite_Dmin.txt"
#infile <- "Dsuite_outs/all-ccgp-nohatch-nobad-trout_dsuite/all-ccgp-nohatch-nobad-trout-with-Chinook-outgroup_dsuite_Dmin.txt"
#infile <- "snphylo-pop-assign_all-ccgp-trout-with-Chinook_dsuite_Dmin.txt"
#infile <- "Dsuite_outs/GL-snphylo-pop-assign_all-ccgp-trout-with-Chinook_dsuite_Dmin.txt"
#infile <- "Dsuite_outs/GL-snphylo-pop-assign_all-ccgp-trout-with-Chinook_dsuite_tree.txt"


D_dmin <- read.table(infile,header=T)

# identify trios with high D
plot(D_dmin$Dstatistic, ylab="D",xlab="trio number")
D_dmin[which(D_dmin$Dstatistic > 0.3),]
# identify those with significant p-value
subset(D_dmin,p.value<0.05)
dim(subset(D_dmin,p.value<0.05))
plot(D_dmin$p.value, ylab="p value",xlab="trio number",ylim=c(0,0.05))

D_dmin$p.adjust <- p.adjust(D_dmin$p.value,method="BH")
dim(subset(D_dmin,p.adjust<0.05)) # didn't change lol

# now look at sig adj p-val (<0.05), high D and high f4 ratio
D_dmin_sorted <- D_dmin %>% arrange(-Dstatistic, f4.ratio)
sorted_outfile <- paste0(infile,".sorted.txt")
write.table(D_dmin_sorted, sorted_outfile, quote=F, sep='\t')



# output 

# snphylo-pop-assign_all-ccgp-trout-with-Chinook_dsuite_Dmin.txt
# based on these outputs here's what I'm thinking
# there's some O clarkii clarkii introgression into 
# O mykiss irideus, gairdnerii, and aguabonita populations
# it's all about 9-10%
# so it could be that clarkii --> common ancestor of mykiss
# maybe less likely that clarkii has independently introgressed with all 3
# should definitely look at mitotypes between subspecies/populations
# let's figure out which parts of the genome experienced introgression!
# choosing the following trios for Dinvestigate
