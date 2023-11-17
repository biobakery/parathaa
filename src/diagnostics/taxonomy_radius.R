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
library(ggridges)

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

for(taxa in taxon_set){
  
  taxa_tips <- tree_tips %>% filter(get({{level}})==taxa)
  taxa_nodes <- taxa_tips$node  
  
  distances <- get_all_pairwise_distances(as.phylo(tree), only_clades = taxa_nodes)

  largest_distance <- max(distances)
  taxon_radi[as.character(taxa)] <- largest_distance
}

write.table(unlist(taxon_radi), file=paste0(outputDir, "/", "taxon_radi", ".tsv"))

densityplot(unlist(taxon_radi), xlab="Taxon Radius")

ggsave(file=paste0(outputDir, "/", "Density_with_zeros.png"), device = "png")
#remove 0s and outliers
taxon_radi_no_zero <- unlist(taxon_radi)[-which(unlist(taxon_radi)==0)]
taxon_radi_no_zero <- taxon_radi_no_zero[-which(taxon_radi_no_zero > 2.5*sd(taxon_radi_no_zero))]

densityplot(taxon_radi_no_zero)
ggsave(file=paste0(outputDir, "/", "Density_with_no_zeros.png"), device = "png")



taxon_set_phyla <- tree_tips$Phylum[match(taxon_set, tree_tips$Species)]

plot_data <- data.frame(Species=names(taxon_radi),
                        Phylum=taxon_set_phyla,
                        Radius=unlist(taxon_radi))

remove_phyla <- which(table(plot_data$Phylum) < 40)
remove_phyla
test <- plot_data %>% filter(!Phylum %in% names(remove_phyla)) %>% 
  ggplot(aes(x=Radius, fill=Phylum, y=Phylum)) + geom_density_ridges() +
  theme_bw(base_size = 16)
test
