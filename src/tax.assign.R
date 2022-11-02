#!/usr/bin/env Rscript
require(docopt)
'Usage:
   tax.assign.R [-j <jplace file> -o <output> -t <tree>]

Options:
   -j jplace file with queries placed into reference tree
   -o output directory [default: output]
   -t reference tree with named internal nodes

 ]' -> doc

#To add as arguments: binom error params

opts <- docopt(doc)

library(logging)

# R logging example 
loginfo("Performing analysis data", logger="")

## This function reads in a jplace object and a reference tree with taxonomy assigned
## to interior nodes and outputs taxonomic assignments (with uncertainty) for the reads
## in the jplace object

library(ggtree)
library(tidytree)
library(treeio)
library(dplyr)

## Inputs
jplaceFile <- opts$j #"/Users/mis696/proj/Phylogeny_Taxonomy/V4_merged_aligned_noend.jplace"
in.treeFile <- opts$t 
outDir <- opts$o ##"out" ###
plotTree <- FALSE

## Read in jplace file, tree
in.jplace <- read.jplace(jplaceFile)
query.names <- unique(in.jplace@placements$name)
load(in.treeFile)
in.tree <- resultData$tax_bestcuts

## Index through query sequences to add taxonomy
results <- c()
for(ind in query.names){
  query.place.data <- in.jplace@placements %>% filter(name==ind)
  tree.w.placements <- left_join(in.tree, query.place.data, by="node")
  tree.w.placements <- tree.w.placements
  tree.w.placements$placementNODE <- NA
  tree.w.placements$placementNODE[tree.w.placements$distal_length>=tree.w.placements$branch.length/2 & !is.na(tree.w.placements$distal_length)] <- 
    tree.w.placements$parent[tree.w.placements$distal_length>=tree.w.placements$branch.length/2 & !is.na(tree.w.placements$distal_length)]
  tree.w.placements$placementNODE[tree.w.placements$distal_length<tree.w.placements$branch.length/2 & !is.na(tree.w.placements$distal_length)] <- 
    tree.w.placements$node[tree.w.placements$distal_length<tree.w.placements$branch.length/2 & !is.na(tree.w.placements$distal_length)]
  
  maxLRs <- tree.w.placements %>% group_by(placementNODE) %>% 
    summarize(maxLR = max(like_weight_ratio, na.rm=T)) %>%
    rename(node=placementNODE)
  tree.w.placements <- left_join(tree.w.placements, maxLRs)
  tree.w.placements$indicator <- NA
  tree.w.placements$indicator[!is.na(tree.w.placements$distal_length)] <- TRUE
  
  tree.w.placements$isPlaced<-NA
  tree.w.placements$isPlaced[tree.w.placements$node %in% tree.w.placements$placementNODE] <- TRUE
  #tree.w.placements$like_weight_ratio[which(tree.w.placements$isPlaced)]
  
  if(plotTree==TRUE){
    mrca <- MRCA(tree.w.placements, query.place.data$node)
    query.subset.tree <- tree_subset(as.treedata(tree.w.placements), mrca$node, levels_back=1)
    query.subset.tree <- as_tibble(query.subset.tree)
    p <- ggtree(as.treedata(query.subset.tree))+ geom_tippoint(aes(color=Genus, shape=isPlaced, size=maxLR)) + 
    geom_nodepoint(aes(color=Genus, shape=isPlaced, size=maxLR))
    ggsave(p, filename = file.path(outDir, ind, "_tree.png"))
  }
  
  assignment <- tree.w.placements %>% filter(maxLR==max(maxLR, na.rm = TRUE))
  assignment$query.name <- ind
  results <- bind_rows(results, assignment)
}

write.table(results, file=file.path(opts$o, "taxonomic_assignments.tsv"),sep = '\t', quote = F, row.names=F)


