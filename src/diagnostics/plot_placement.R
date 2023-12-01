## ids: names of query sequences
## in.tree.data: reference tree with interior nodes named (treedata)
## in.jplace: jplace object with placements
## compare.synth: from synthetic.benchmark.R. needs a column called "Genus.silva" and "Genus.parathaa" (or whatever level you use)
##    to extract the true and assigned names from, and "AccID" to match the query ID
## level: level to examine and plot taxonomy
## outputDir: Output directory name

##First, get unique interior node labels:
# in.tree.data <- in.tree.data %>% mutate(label_new = ifelse(isTip, label, paste0("Node_",node))) %>%
#   select(!label) %>%
#   rename(label=label_new)



plot_placement <- function(ids, in.tree.data, in.jplace, compare.synth, level="Genus", outputDir="output", levels_back=3){

  for(nm in ids){ 
    plotTree <- as.phylo(in.tree.data)
    plc <- in.jplace@placements[which(in.jplace@placements$name==nm),]
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
      truName <- compare.synth %>% filter(AccID==nm) %>% select(paste0(level,".silva")) %>% as.character()
      parName <- compare.synth %>% filter(AccID==nm) %>% select(paste0(level,".parathaa")) %>% as.character()
      
      
      levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      plevel <- level
      if(length(table(plotTree2[[level]]))>12)
        plevel <- levels[which(levels=="Genus")-1]
      plots.out[[paste0(nm, "_placement", pind)]] <-  ggtree(as.treedata(plotTree2), aes_string(color=plevel)) + geom_tippoint() + geom_nodepoint() + geom_tiplab() + 
        labs(title=paste0("Source: ", truName),  subtitle = paste0("Parathaa: ", parName)) +
        scale_size_manual(values=c(1, 5)) + geom_treescale()
      ggsave(paste0(outputDir, "/", nm, "_placement", pind, ".png"), 
             plot=plots.out[[paste0(nm, "_placement", pind)]] ,
             width = 12, height=12, units = "in")
      
    }
  }
  return(plots.out)
}

# Example:
##plot_placement(ids=compare.synth$AccID[1:5], in.tree.data=in.tree.data, in.jplace=in.jplace, compare.synth=compare.synth, outputDir="test_function")