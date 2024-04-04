## trout subspecies tree
library(ape)
library(ggtree)
library(ggplot2)

# set trout palette
trout_palette <- c("#96ACC0", "#88716C", "#F6C7C0", "#324141", "#DDD4DF","#BFD8D4","#6F6B4B")

# read metadata
metadata <- read.table("refs/SETS-w-Chinook.plus.tsv",header=T)

# Load in the newick file
ml_tree <- read.tree("snphylo_dnaml_trees/run2_l0.2M0.4/name-swapped.snphylo.l0.2.M0.4.out.bs.tree")
outgroup <- subset(metadata,species_name=="Outgroup")$new_sample_name
ml_tree <- root(ml_tree, outgroup) 


tt <- ggtree(ml_tree) + geom_tiplab(size=2) + ggplot2::xlim(0,0.5) +
  geom_nodelab(aes(label=label, subset=as.numeric(label) > 75),size=2,hjust=-0.01)
#  geom_strip('gairdnerii', 'aguabonita', barsize=2, color=trout_palette[1], 
#             label="O. mykiss", offset=1.5, offset.text=.2) +
#  geom_strip('seleniris', 'clarkii', barsize=2, color=trout_palette[2], 
#             label="O. clarkii", offset=1.5, offset.text=.2) 
# %<+% metadata +
#  geom_tippoint(aes(color=Country)) + # color code tips
#  geom_treescale(x=0, y=45, fontsize=4, linesize=2, offset=2, width=10) # add scale

plot(tt)
ggsave(tt,filename="Figures/chinook-rooted.snphylo.l0.2.M0.4.out.bs.pdf",height=10,width=10,bg="transparent")


# color tips
x <- full_join(as_tibble(ml_tree), metadata, by = join_by("label" == "new_sample_name"))
meta_matched_ml_tree <- as.treedata(x)

trout_palette <- c("#96ACC0", "#88716C", "#F6C7C0", "#324141", "#DDD4DF","#BFD8D4","#6F6B4B","#8DA750","salmon","navy")
tt_color <-
  ggtree(meta_matched_ml_tree) + ggplot2::xlim(0,0.5) +
  geom_nodelab(aes(label=label, subset=as.numeric(label) > 75),size=2,hjust=-0.01) +
  geom_tiplab(aes(color = subspecies_name),size=2) +
  scale_color_manual(values=trout_palette)
ggsave(tt_color,filename="Figures/colorful-tree.chinook-rooted.snphylo.l0.2.M0.4.out.bs.pdf",height=10,width=10,bg="transparent")


