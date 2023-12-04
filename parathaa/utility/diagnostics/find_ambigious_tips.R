### This script finds sequences within the reference database that if classified by Parathaa using its 
### optimal threshold it would be ambigious at that taxonomic level
### This gives an idea to the user of what taxa are amibigious within the region that the user choose to use

require(docopt)

'Usage:
  find_ambigious_tips.R [-t <annotated.tree> -l <taxa_level> -o <output_dir> -d <distance> -s <TRUE_FALSE]
  
  
Options:
  -t annotated tree file created by Parathaa run_tree_analysis.py
  -o output directory
  -l taxonomic level to check for disagreements at
  -d allowed distance between tips to call them ambigious
  -s summarize TRUE/FALSE
  ]' -> doc


library(dplyr)
library(castor)
library(treeio)

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

threshold <- opts$d

summarize <- opts$s

Ambig_tips <- list()
#go through each tip

for(node in tree_tips$node){
  distances <- get_pairwise_distances(as.phylo(tree), rep(node, length(tree_tips$node)), tree_tips$node)

  amibg_nodes <- which(distances <= threshold)
  
  if(length(amibg_nodes) <= 1){
     next
  }
  
  ambig_taxa_list <- unlist(tree_tips[amibg_nodes,level])
  
  if(length(unique(unlist(ambig_taxa_list))) <= 1){
    next
  }
    
  Ambig_tips[[as.character(node)]] <- ambig_taxa_list
  
}

seq_names <- tree_tips$label[match(names(Ambig_tips),tree_tips$node)]

names(Ambig_tips) <- seq_names

if(length(seq_names)==0){
  message("No ambigious tips")
  quit()
}

max.length <- max(sapply(Ambig_tips, length))
Ambig_tips <- lapply(Ambig_tips, function(v) { c(v, rep(NA, max.length-length(v)))})
multi_match_df <- data.frame(do.call(rbind, Ambig_tips))

colnames(multi_match_df) <- paste("Taxa", seq(1, length(colnames(multi_match_df))), sep="")
multi_match_df[is.na(multi_match_df)] <- ""

write.table(multi_match_df, file=paste0(outputDir, "/", "ambigious_tips_d", threshold, ".tsv"), sep = "\t",
            col.names = NA, row.names = TRUE)



ambig_taxa <- tree_tips[match(rownames(multi_match_df), tree_tips$label),level]


#write summary stats (i.e precent of sequences that are from that species that are ambigious)
if(summarize){
  n_seqs_ambig <- table(ambig_taxa)
  summary_df <- data.frame(Taxa=names(n_seqs_ambig),
                           N_seqs=as.vector(n_seqs_ambig))

  n_taxa <- table(tree_tips[,level])  
  percent <- n_seqs_ambig/n_taxa[match(names(n_seqs_ambig), names(n_taxa))]
  
  summary_df$percent <- percent
  
  write.table(summary_df, file=paste0(outputDir, "/", "ambigious_tips_summary_d", threshold, ".tsv"), sep = "\t",
              col.names = NA, row.names = TRUE, quote=FALSE)
}




