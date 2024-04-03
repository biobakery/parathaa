require(docopt)

'Usage:
  plot_assignment.R [--parathaa_db_tree <parathaa_db> --level <taxonomic level> --steps_back <integer> --assignments <parathaa_assignments> --jplace <jplace_file> -o <output_dir> --id <query_id>]
  
  
Options:
  --parathaa_db_tree parathaa database
  --assignments parathaa assignments
  --jplace jplace file
  --level taxonomic level that tree is labelled with
  --steps_back the number of steps back from the placement tip that the root should be
  --id query id to plot
  -o output
  ]' -> doc

opts <- docopt(doc)

library(treeio)
library(dplyr)
library(phytools)
library(ggtree)
library(ggplot2)

if(!interactive()) pdf(NULL)

load(opts$parathaa_db_tree)
in.tree <- resultData$tax_bestcuts

in.tree <- in.tree %>% mutate(label_new = ifelse(isTip, label, paste0("Node_",node))) %>%
  select(!label) %>%
  rename(label=label_new)

jplaceFile <-  opts$j
in.jplace <- read.jplace(jplaceFile)

query <- opts$id

assignments <- read.table(opts$assignments, sep="\t", header=T, row.names=1)


plot_placement <- function(ids, in.tree.data, in.jplace, level="Genus", outputDir="output", levels_back=3, assignments){
  
  for(nm in ids){ 
    plotTree <- as.phylo(in.tree.data)
    plc <- in.jplace@placements[which(in.jplace@placements$name==nm),]
    #filter to only those that are used for assignment
    plc <- plc %>% filter(like_weight_ratio >= max(plc$like_weight_ratio)*0.5)
    plots.out <- list()
    for(pind in 1:nrow(plc)){
      ## Add placement to reference tree
      plotTree <- as.phylo(in.tree.data)
      plotTree <- bind.tip(plotTree, tip.label=paste0("Placement_", pind), edge.length=plc$pendant_length[pind], where=plc$node[pind], position=plc$distal_length[pind])
      plc1 <- which(plotTree$tip.label==paste0("Placement_", pind))
      ## Subset tree for plotting
      test.tree <- tree_subset(plotTree, node=plc1, levels_back = levels_back)
      ## For plotability, only go 2 levels back if >50 tips
      if(nrow(as_tibble(test.tree))>50)
        test.tree <- tree_subset(plotTree, node=plc1, levels_back = levels_back-1)
      
      ## Re-merge tree with taxonomy info from reference tree by label
      labelInfo <- in.tree.data %>% select(label, Kingdom, Phylum, Class, Order, Family, Genus, Species)
      plotTree2 <- left_join(as_tibble(test.tree), labelInfo)
      ## Get real taxonomy and assigned taxonomy for query
      parName <- assignments[nm,-8]
      Name <- paste(unlist(parName), collapse=";")
      
      levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      plevel <- level
      if(length(table(plotTree2[[level]]))>12)
        plevel <- levels[which(levels=="Genus")-1]
      plots.out[[paste0(nm, "_placement", pind)]] <-  ggtree(as.treedata(plotTree2), aes_string(color=plevel)) + geom_tippoint() + geom_nodepoint() + geom_tiplab() + 
        labs(title=paste0("Assignment: ", Name)) +
        scale_size_manual(values=c(1, 5)) + geom_treescale()
      ggsave(paste0(outputDir, "/", nm, "_placement", pind, ".png"), 
             plot=plots.out[[paste0(nm, "_placement", pind)]] ,
             width = 12, height=12, units = "in")
      
    }
  }
  return(plots.out)
}

plot_placement(ids = query, in.tree.data = in.tree, in.jplace = in.jplace, level = opts$level, outputDir = opts$o, levels_back = opts$steps_back, assignments = assignments)
