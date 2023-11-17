### This script takes in an annotated tree and then outputs tips that have a taxonomic label 
### that do not agree with
### its parent node at the given level specified. 
### This script was used to QC the MetaRef 16S database to remove potential cases where mis-assembly or
### mis-binning resulted in 16S sequences in isolates that were not actually from that isolate


require(docopt)

'Usage:
  Id.tip.disagreement.R [-t <annotated.tree> -l <taxa_level> -p <TRUE_FALSE> -o <output_dir>]
  
  
Options:
  -t annotated tree file created by Parathaa run_tree_analysis.py
  -o output directory
  -l taxonomic level to check for disagreements at
  -p Whether to plot the subregions of tip disagreements
  ]' -> doc




opts <- docopt(doc)

library(ggtree)
library(ggplot2)
library(dplyr)
library(treeio)

#plotting function
generate_subtree <- function(tree, node, levels_back=3, level, auto_set=TRUE){
  subtree <- tree_subset(as.phylo(tree), node=node, levels_back = levels_back)
  subtree_tibble <- as_tibble(subtree)
  subtree_tibble <- subtree_tibble %>% left_join(tree[,c(10:15, 4)], by='label')
  
  if(auto_set){
    n_tips_okay <- FALSE
  }else{
    n_tips_okay <- TRUE
  }
  
  #move back a level if more than 50 tips or if the taxonomic level has too many uniques
  while(!n_tips_okay){
    #message("Setting levels back auto")
    if(levels_back==1){
      n_tips_okay <- TRUE
    }else if(nrow(as_tibble(subtree))>50 | length(unique(subtree_tibble[,level])) > 12){
      #message(levels_back)
      levels_back <- levels_back-1
      subtree <- tree_subset(as.phylo(tree), node=node, levels_back = levels_back)
      subtree_tibble <- as_tibble(subtree)
      subtree_tibble <- subtree_tibble %>% left_join(tree[,c(10:15, 4)], by='label')
    }else{
      n_tips_okay <- TRUE
    }
  }
  
  
  #message(levels_back)
  level <- sym(level)
  
  plot <- ggtree(as.phylo(subtree_tibble)) %<+% subtree_tibble + geom_tiplab(aes(color=!!level)) + 
    geom_nodepoint(aes(color=!!level)) +
    ggtitle(paste(tree$label[which(tree$node==node)])) + geom_treescale()
  
  return(plot)
}


load(opts$t)

tree <- resultData$tax_bestcuts

level <- as.character(opts$l)

level_parent <- paste("parent_",level, sep="")

message(level)

plot <- opts$p
outputDir <- opts$o

## Prep tree
tree$label <- make.unique(tree$label)
tree <- tree %>% arrange(node)


tree_tips <- tree[which(tree$isTip),]



tree_tips[[level_parent]] <- unlist(tree[tree_tips$parent,level])


agreement_index <- apply(tree_tips, 1, function(x) grepl(x[level], x[level_parent]))




tip_disagreement <- tree_tips[which(agreement_index==FALSE),]

#remove NA parents
tip_disagreement <- tip_disagreement[!is.na(tip_disagreement[,level_parent]),]
tip_disagreement <- tip_disagreement[!tip_disagreement[,level_parent]=="",]


if(plot){
  plot_list <- list()
  for(node in tip_disagreement$node){
    #message(node)
    plot_list[[as.character(node)]] <- generate_subtree(tree, node, levels_back = 3, level, auto_set = TRUE)
  }
  for(n in plot_list){
    ggsave(paste0(outputDir, "/", n$labels$title, ".png"), plot=n, width=12, height=12, unit="in")
  }

}

tip_list <- write.table(tip_disagreement$label, file = paste0(outputDir, "/", "non_agreement_tips.tsv"),
                        col.names=FALSE, row.names = FALSE, sep="\t", quote=FALSE)

write.table(tip_disagreement, file=paste0(outputDir, "/", "non_agreement_fill_info.tsv"), sep = "\t",
            col.names = NA, row.names = TRUE)


