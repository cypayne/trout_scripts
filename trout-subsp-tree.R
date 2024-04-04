## trout subspecies tree
library(ape)
library(ggtree)
library(ggplot2)

# set trout palette
trout_palette <- c("#96ACC0", "#88716C", "#F6C7C0", "#324141", "#DDD4DF","#BFD8D4","#6F6B4B","#8DA750","salmon","navy")


## simple subspecies tree
# Load in the newick file
#trout_tree <- read.tree("trout-subspecies-tree.ml.tre")
trout_tree <- read.tree("snphylo_dnaml_trees/run7_l0.2.M0.25.p25_no-hatchery-no-bad/snphylo.l0.2.M0.25.p25_broad-group.ml.tre")

trout_tree <- root(trout_tree, 'Outgroup') 

tt <- ggtree(trout_tree) + geom_tiplab(size=5) + ggplot2::xlim(0,10) +
      geom_strip('gairdnerii', 'aguabonita', barsize=2, color=trout_palette[1], 
           label="O. mykiss", offset=1.5, offset.text=.2) +
      geom_strip('seleniris', 'clarkii', barsize=2, color=trout_palette[8], 
             label="O. clarkii", offset=1.5, offset.text=.2) 
tt

#ggsave(tt,filename="Figures/simplified_ccgp-trout_subspecies_tree.pdf",height=10,width=10,bg="transparent")


d <- data.frame(node=c(10, 19), species=c("mykiss", "clarkii"))
tt <- ggtree(trout_tree) + geom_tiplab(size=5) + ggplot2::xlim(0,10) + 
      geom_hilight(data=d, aes(node=node, fill=species),
               type = "roundrect")
tt



## 

library(ggrepel)
library("ggtree")
library("treeio")

metadata <- read.csv("/Users/cypayne/Desktop/trout_projects/refs/CCGP-trout_sample-coordinates_correctly-geo-ranked_with-Chinook-Outgroup.csv",header=T)
tree_file <- "/Users/cypayne/Desktop/trout_projects/snphylo_dnaml_trees/run7_l0.2.M0.25.p25_no-hatchery-no-bad/rooted.name-swapped.snphylo_no-hatchery_no-bad.l0.2.M0.25.p25.out.bs.tre"
#tree <- read.tree("snphylo_dnaml_trees/run7_l0.2.M0.25.p25_no-hatchery-no-bad/name-swapped.snphylo_no-hatchery_no-bad.l0.2.M0.25.p25.out.bs.tree")
tree <- read.tree(tree_file)
tree_labels <- as.data.frame(tree$tip.label)
colnames(tree_labels) <- c("dsuite_sample_name")
tree_metadata <- merge(tree_labels,metadata,by="dsuite_sample_name")
# replace tree tip labels with new name
tree$tip.label <- tree_metadata$coord_name

ggplot(tree, aes(x, y)) + geom_tree(layout="fan") + theme_tree()
ggtree(tree, layout="fan", open.angle=120) 

ggtree(tree) + geom_label_repel(aes(label=bootstrap, fill=bootstrap)) + 
  theme(legend.position = c(.1, .8)) + scale_fill_viridis_c()

ggtree(tree, layout="fan", open.angle=90) + geom_label2(aes(label=label, fill=as.numeric(label),
                             subset = !is.na(as.numeric(label)) & as.numeric(label) > 80), label.size=0.15) +
  theme(legend.position = c(.1, .8)) + scale_fill_viridis_c()

tree_plot <-
  ggtree(tree, layout="fan", open.angle=150) + 
  geom_nodepoint(aes(fill=as.numeric(label), subset = !is.na(as.numeric(label)) & as.numeric(label) > 80), size=2) +
  theme(legend.position = c(.1, .8)) + scale_fill_continuous(low="blue",high="yellow") +
  geom_tiplab() 
tree_plot



geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint(aes(shape = trophic_habit, color = trophic_habit, 
                    size = mass_in_kg)) + 
  theme(legend.position = "right") + 
  scale_size_continuous(range = c(3, 10))

p2 %<+% df_inode_data + 
  geom_label(aes(label = vernacularName.y, fill = posterior)) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu"))
