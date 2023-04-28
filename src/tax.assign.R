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
library(phytools)

## Inputs
jplaceFile <-  opts$j
in.treeFile <-  opts$t
outDir <- opts$o 
plotTree <- FALSE

## Read in jplace file, tree
in.jplace <- read.jplace(jplaceFile)
query.names <- unique(in.jplace@placements$name)
load(in.treeFile)
in.tree <- resultData$tax_bestcuts

## Index through query sequences to add taxonomy
results <- c()
# might beable to speed this up using foreach/ vectorized code
for(ind in query.names){
  ## Read in data, filter to most likely placement(s) and make various useful formats of it
  query.place.data <- in.jplace@placements %>% filter(name==ind)  %>% filter(like_weight_ratio==max(like_weight_ratio))
  tree.w.placements <- left_join(in.tree, query.place.data, by="node")
  tree.w.placements.phy <- tidytree::as.phylo(tree.w.placements)
  tree.w.placements.tib <- as_tibble(tree.w.placements)
  
  ## Add a column of heights for each node
  tree.w.placements.tib$nodeHeight <- 0
  
  # we want to see how close that sequence is to the children node 
  # we could have something that all "original" children as X genus but then we place new node very far away which deosn't mean its that genus
  tree.w.placements.tib$nodeHeight[tree.w.placements.phy$edge[,"node"]] <- nodeHeights(tree.w.placements.phy)[,2]
  
  ## Identify query placements and their offspring, extract maximum node heights
  ind.offs <- offspring(tree.w.placements.phy, query.place.data$node, self_include=T) ## query.place.data$node might be a vector... 
  
  
  # check for tied placements
  if(length(query.place.data$node)==1){
    maxNodeHeights <- tree.w.placements.tib %>% filter(node == query.place.data$node) %>% summarize(max(nodeHeight)) %>% as.numeric
  } else {
  maxNodeHeights <- sapply(ind.offs, FUN= function(x) tree.w.placements.tib %>% filter(node %in% x) %>% summarize(max(nodeHeight)) %>% as.numeric)
  }
  
  # calculate the distance from the query node to any other child node
  maxDistPlacements <- maxNodeHeights-tree.w.placements.tib$nodeHeight[query.place.data$node] + tree.w.placements.tib$distal_length[query.place.data$node] + tree.w.placements.tib$pendant_length[query.place.data$node]
  
  if(length(maxDistPlacements)==1) ## If only one placement, maxDistPlacements isn't named, so add name
    names(maxDistPlacements) <- query.place.data$node
  
  ## Load thresholds for levels
  load(file.path(opts$o, "optimal_scores.RData"))
  bestThresh <- plotData2 %>% group_by(Level) %>% summarise(minThreshold = mean(minThreshold))
  cutoffs <- bestThresh$minThreshold
  names(cutoffs) <- bestThresh$Level
  
  ## Find lowest justifiable taxonomic classification of placement 
  ## Using max of distances across placements to be conservative, for now
  assignmentLevels <- names(which(cutoffs>max(maxDistPlacements)))
  assignment <- tree.w.placements.tib[query.place.data$node, c("Kingdom", rev(assignmentLevels))]
  ## Deal with multifurcation or placements at child nodes
  #this is the case when multiple tides are below a single node?
  #or is this dealing with a case where p placer places a sequence as a child of a placed sequence?
  if(length(unique(tree.w.placements.tib$parent[query.place.data$node]))==1 & length(query.place.data$node)>1)
    assignment <- tree.w.placements.tib[unique(tree.w.placements.tib$parent[query.place.data$node]), c("Kingdom", rev(assignmentLevels))]

  
  assignment$query.name <- ind
  results <- bind_rows(results, assignment)
}

pick.taxon <- function(x){
  if(sum(!is.na(x))==0){
    x3 <- NA
  } else{
    x2 <- sort(unique(unlist(strsplit(x, ";"))))
    x3 <- paste(x2, collapse = ";")
  }
  return(x3)
}

tax_parathaa <- results %>%
  group_by(query.name) %>% 
  summarize(Kingdom  = pick.taxon(Kingdom),
            Phylum = pick.taxon(Phylum),
            Class = pick.taxon(Class),
            Order = pick.taxon(Order),
            Family = pick.taxon(Family),
            Genus = pick.taxon(Genus),
            Species = pick.taxon(Species))

write.table(tax_parathaa, file=file.path(opts$o, "taxonomic_assignments.tsv"),sep = '\t', quote = F, row.names=F)


