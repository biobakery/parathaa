#!/usr/bin/env Rscript
require(docopt)
'Usage:
   tax.assign.R [-j <jplace file> -o <output> -t <tree> -s <optimal_scores> --threads <threads> -d <delta>]

Options:
   -j jplace file with queries placed into reference tree
   -o output directory [default: output]
   -t reference tree with named internal nodes
   -s Optimal threshold scores
   --threads number of threads to run in parallel [default: 1]
   -d delta [default: 0.02]

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
library(doSNOW)
source("utility/nearest.neighbor.revisions.R")


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Inputs
jplaceFile <-  opts$j
in.treeFile <-  opts$t
optimalFile <- opts$s
outDir <- opts$o 
plotTree <- FALSE

## Read in jplace file, tree
in.jplace <- read.jplace(jplaceFile)
query.names <- unique(in.jplace@placements$name)
load(in.treeFile)
load(optimalFile)
in.tree <- resultData$tax_bestcuts

in.tree <- in.tree %>% mutate(label_new = ifelse(isTip, label, paste0("Node_",node))) %>%
  select(!label) %>%
  rename(label=label_new)

delta <- as.numeric(opts$d)


bestThresh <- plotData2 %>% group_by(Level) %>% summarise(minThreshold = mean(minThreshold))
cutoffs <- bestThresh$minThreshold
names(cutoffs) <- bestThresh$Level

## Index through query sequences to add taxonomy
# might beable to speed this up using foreach/ vectorized code
iterations <- length(query.names)
pb <- txtProgressBar(max=iterations, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
pb_opts <- list(progress=progress)

cl <- makeCluster(as.numeric(opts$threads))
registerDoSNOW(cl)

result <- foreach(i=1:length(query.names), .combine = bind_rows, .options.snow = pb_opts,
                  .packages = c("dplyr", "treeio", "tidytree", "phytools")) %dopar% {
  
  ind <- query.names[i]
  setTxtProgressBar(pb,i)
  ## Read in data, filter to most likely placement(s) and make various useful formats of it
  query.place.data <- in.jplace@placements %>% filter(name==ind)  %>% filter(like_weight_ratio > 0.5*max(like_weight_ratio))
  tree.w.placements <- left_join(in.tree, query.place.data, by="node")
  class(tree.w.placements) <- c("tbl_tree", class(tree.w.placements))
  tree.w.placements.phy <- tidytree::as.phylo(tree.w.placements)
  tree.w.placements.tib <- as_tibble(tree.w.placements)
  
  ## Add a column of heights for each node
  tree.w.placements.tib$nodeHeight <- 0
  
  # we want to see how close that sequence is to the children node 
  # compute the height above the root for all nodes
  tree.w.placements.tib$nodeHeight[tree.w.placements.phy$edge[,"node"]] <- nodeHeights(tree.w.placements.phy)[,2]
  
  ## Identify query placements and their offspring, extract maximum node heights
  ind.offs <- offspring(tree.w.placements.phy, query.place.data$node, self_include=T) ## query.place.data$node might be a vector... 
  
  
  # check for tied placements
  # check for tied placements
  if(length(query.place.data$node)==1){
    ##maxNodeHeights <- tree.w.placements.tib %>% filter(node == query.place.data$node) %>% summarize(max(nodeHeight)) %>% as.numeric ## old
    maxNodeHeights <- tree.w.placements.tib %>% filter(node %in% ind.offs) %>% summarize(max(nodeHeight)) %>% as.numeric ## new
  } else {
    maxNodeHeights <- sapply(ind.offs, FUN= function(x) tree.w.placements.tib %>% filter(node %in% x) %>% summarize(max(nodeHeight)) %>% as.numeric)
  }
  
  # calculate the distance from the query node to any other child node
  #maxNodeHeights is distance from the root to the tip
  ## shouldn't this really be the distance from the furthest child node and the root? 
  ## I don't think this is being calculated correctly when sequences are placed on internal nodes which is messy up
  ## the assignment algorithim.
  #tree.w.placements.tin$nodeHeight[query.place.data$node]
  ## this is the distance from the root to the query node
  
  #it seems like maxNodeHeights
  maxDistPlacements <- maxNodeHeights-tree.w.placements.tib$nodeHeight[query.place.data$node] + tree.w.placements.tib$distal_length[query.place.data$node] + tree.w.placements.tib$pendant_length[query.place.data$node]
  
  #what does this do?
  if(length(maxDistPlacements)==1) ## If only one placement, maxDistPlacements isn't named, so add name
    names(maxDistPlacements) <- query.place.data$node
  
  ## Load thresholds for levels

  
  ## Find lowest justifiable taxonomic classification of placement 
  ## Using max of distances across placements to be conservative, for now
  numLevels <- lapply(maxDistPlacements, FUN=function(x) names(which(cutoffs>x))) %>% lapply(length) %>% getmode %>% unlist
  assignmentLevels <- names(cutoffs[1:numLevels]) ### LEFT OFF HERE

  assignment <- tree.w.placements.tib[query.place.data$node, c("Kingdom", assignmentLevels)]

  
  
  #if there is a multifurcation then assign the taxonomy as its parent node
  if(length(unique(tree.w.placements.tib$parent[query.place.data$node]))==1 & length(query.place.data$node)>1)
    assignment <- tree.w.placements.tib[unique(tree.w.placements.tib$parent[query.place.data$node]), c("Kingdom", rev(assignmentLevels))]

  
  assignment$query.name <- ind
  assignment$maxDist <- max(maxDistPlacements)
  #assignment$same_testing <- maxNodeHeights == tree.w.placements.tib$nodeHeight[query.place.data$node]
  #I think it both cases the calculations are only correct for tips and not internal node placements.
  return(assignment)
}
## what is the point of this function??
## oh is this for when there is multiple assignments?

pick.taxon <- function(x){
  if(sum(!is.na(x))==0){
    x3 <- NA
  } else{
    x2 <- sort(unique(unlist(strsplit(x, ";"))))
    x3 <- paste(x2, collapse = ";")
  }
  return(x3)
}

tax_parathaa <- result %>%
  group_by(query.name) %>% 
  summarize(Kingdom  = pick.taxon(Kingdom),
            Phylum = pick.taxon(Phylum),
            Class = pick.taxon(Class),
            Order = pick.taxon(Order),
            Family = pick.taxon(Family),
            Genus = pick.taxon(Genus),
            Species = pick.taxon(Species))


## nearest.neighbor.distance code is extremely slow need to fix this at somepoint
## not a prority issue though.
if(delta>0){
  dists <- nearest.neighbor.distances(tax.df=tax_parathaa, 
                                      placement.object=in.jplace, 
                                      reference.tree=in.tree, 
                                      max.radius=delta)

  tax_parathaa <- nearest.neighbor.revisions(tax.df=tax_parathaa, 
                                             distances=dists, 
                                             radius=delta)
}
tax_parathaa$maxDist <- result$maxDist[match(tax_parathaa$query.name, result$query.name)]

write.table(tax_parathaa, file=file.path(opts$o, "taxonomic_assignments.tsv"),sep = '\t', quote = F, row.names=F)

