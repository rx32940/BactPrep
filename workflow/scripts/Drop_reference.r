library(ggtree)
library(dplyr)
library(treeio)

# read Gubbins Final tree
tempest_tree <- read.newick("/Users/rx32940/Downloads/Lint123467.final_tree.tre")

# plot tree
tempest_tree_p <- ggtree(tempest_tree)+
  geom_text(aes(label=node,size=12))+
  geom_tiplab()

# get the tree tip "Reference"
tip_2_drop <- tempest_tree_p$data %>% subset(label %in% c( "Reference" ) )

# Drop reference tip 
subset1_drop <- drop.tip(tempest_tree, tip_2_drop$node)

# read reference dropped tree into a new tree
subset1_tree <- ggtree(subset1_drop)+
  geom_text(aes(label=node))+
  geom_tiplab()


# write the tree
write.tree(as.phylo(subset1_tree),"/Users/rx32940/Downloads/Lint123467.final_tree.L23DropOutlierDropped.tre")