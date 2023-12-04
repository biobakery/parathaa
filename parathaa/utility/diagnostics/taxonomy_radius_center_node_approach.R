#Script that will caluculate the species radius for each species in the phyologentic tree


require(docopt)

'Usage:
  taxonomy_radius.R [-t <annotated.tree> -l <taxa_level> -o <output_dir>]
  
  
Options:
  -t annotated tree file created by Parathaa run_tree_analysis.py
  -o output directory
  -l taxonomic level to check for disagreements at
  ]' -> doc


library(castor)
library(dplyr)
library(treeio)
library(ggplot2)
library(ggtree)

opts <- docopt(doc)

load(opts$t)

tree <- resultData$tax_bestcuts

level <- as.character(opts$l)
message(level)
outputDir <- opts$o

## Prep tree
tree$label <- make.unique(tree$label)
tree <- tree %>% arrange(node)

tree_tips <- tree[which(tree$isTip),]

taxon_set <- unique(unlist(tree_tips[,level]))
taxon_set <- taxon_set[-which(is.na(taxon_set))]

taxon_radi <- list()
tree$isCenter <- FALSE
center_nodes <- c()
for(taxa in taxon_set){
  
  taxa_tips <- tree_tips %>% filter(get({{level}})==taxa)
  taxa_nodes <- taxa_tips$node  
  
  if(length(taxa_nodes) == 1){
    taxon_radi[as.character(taxa)] <- 0
  }
  
  distances <- get_all_pairwise_distances(as.phylo(tree), only_clades = taxa_nodes)
  test <- apply(distances, 1, max, na.rm=TRUE)
  center_tip <- taxa_nodes[which.min(test)]
  center_nodes <- c(center_nodes, center_tip)
  tree$isCenter[which(tree$node==center_tip)] <- TRUE
  #node 7620 is the center node
  #next we extract a subtree with maximum radius around the center node
  #okay we have the extracted tree radius now we need to find the optimal distance
  #whats the best way of choosing the maximum distance range...
  #there is probably a better way to determine the distances to try here...
  cutoffs<-c(seq(0.001, 0.01, by=0.001), seq(0, max(distances), by=0.1*max(distances)))
  cutoffs <- sort(cutoffs)
  cut_res <- c()
  for(cut in cutoffs){
    tree_radius <- extract_tip_radius(as.phylo(tree), focal_tip = center_tip, radius = cut, include_subtree = F)
    
    #inclusion score is number of tips that are 
    taxa_tab <- table(tree_tips[match(tree_radius$radius_tips, tree_tips$node),level])
    
    #we want to maximize the inclusion score
    inclusion_score <- taxa_tab[which(names(taxa_tab)==taxa)]/length(taxa_nodes)
    
    #we want to make the exclusion score as low as possible...
    exclusion_score <- sum(taxa_tab[which(names(taxa_tab)!=taxa)])/length(taxa_nodes)
    
    total_score <- inclusion_score - exclusion_score
    cut_res[as.character(cut)] <- total_score
  }
  
  

  largest_distance <- as.numeric(names(cut_res)[which.max(cut_res)])
  taxon_radi[as.character(taxa)] <- largest_distance
}

#write raw table of optimal radi
write.table(unlist(taxon_radi), file=paste0(outputDir, "/", "taxon_radi", ".tsv"))

tree$optimal_radi <- NA
tree$optimal_radi[center_nodes] <- unlist(taxon_radi)
## plot across the tree?


plot_df <-  tree %>% select(label, optimal_radi, Phylum, isCenter)
phyla_n <- plot_df %>% group_by(Phylum) %>% summarize(n())

Other_phlya <- phyla_n$Phylum[which(phyla_n$`n()` < 50)]

plot_df$Phyla_agg <- plot_df$Phylum
plot_df$Phyla_agg[which(plot_df$Phyla_agg %in% Other_phlya)] <- "Other"
plot_df$test2 <- plot_df$optimal_radi
plot_df$phyla2 <- plot_df$Phyla_agg
plot_df2 <- data.frame(plot_df)
plot_df2$Phyla <- plot_df2$Phyla_agg

tree_plot <- ggtree(as.phylo(tree))

full_plot <- tree_plot %<+% tree + geom_tippoint(aes(alpha=isCenter), size=3) + 
  scale_alpha_manual(values=c(0,0.5))


facet_plot(full_plot, 
           panel="Log Optimal Radi", data=plot_df2, mapping=aes(x=log(test2+0.01), color=Phyla), geom=geom_point) +
  theme_bw(base_size=10)


ggsave(filename = paste0(outputDir, "/", "tree_with_taxon_radi.png"), device = "png", width=12, height=12)

Smooth_plot <- full_plot$data %>% filter(!is.na(x), !is.na(optimal_radi)) %>% 
  ggplot(aes(x=y, y=optimal_radi)) + geom_smooth(color="red") + 
  geom_point(alpha=0.5) + theme_bw(base_size=12) + coord_flip() + xlab("Tree Position") + ylab("Optimal Radii")

Smooth_plot + xlab("Position in tree") + ylab("Optimal Radi")



ggsave(filename = paste0(outputDir, "/", "Optimal_radii_by_tree_position.png"), device="png", width=12, height=12)


Smooth_plot <- full_plot$data %>% filter(!is.na(x), !is.na(optimal_radi), optimal_radi!=0) %>% 
  ggplot(aes(x=y, y=log(optimal_radi))) + geom_smooth() + 
  geom_point(alpha=0.5) + theme_bw(base_size = 12) + coord_flip() + xlab("Tree Position") + ylab("Optimal Radii")

Smooth_plot

ggsave(filename=paste0(outputDir, "/", "Optimal_radii_by_tree_position_no_zero_log.png"), device="png", width=12, height=12)

