#!/usr/bin/env Rscript
require(docopt)
'Usage:
   tax.assign.R [-j <jplace file> -o <output> -t <tree> -s <optimal_scores> --threads <threads> -d <delta> -m <mult> --md <mult_delta> --util1 <nearest_neighbor_PATH> --genus_mult <genus_multi> ]

Options:
   -j jplace file with queries placed into reference tree
   -o output directory [default: output]
   -t reference tree with named internal nodes
   -s Optimal threshold scores
   --threads number of threads to run in parallel [default: 1]
   -d delta [default: 0]
   -m mult [default: 0.1]
   --md [default: 0.5]
   --util1 PATH to nearest_neighbours_parallel.R [default: utility/nearest_neighbours_parallel.R]
   --genus_mult [default: 1]

 ]' -> doc
opts <- docopt(doc)


library(ggtree)
library(tidytree)
library(treeio)
library(dplyr)
library(phytools)
library(doSNOW)
source(opts$util1)


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



### Add species multiplier default is 0.1 in "sensitive mode"
cutoffs["Species"] <- cutoffs["Species"] * as.numeric(opts$m)

### Add genus threshold multipler default is 1 in all modes
cutoffs["Genus"] <- cutoffs["Genus"]*as.numeric(opts$genus_mult)

#if delta multiplier is not left to default (0) change the delta nultipler value (default 0.5)
if(delta!=0){
   delta <- cutoffs["Species"]* as.numeric(opts$md)
}

## Index through query sequences to add taxonomy
# might beable to speed this up using foreach/ vectorized code
iterations <- length(query.names)
pb <- txtProgressBar(max=iterations, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
pb_opts <- list(progress=progress)

cl <- makeCluster(as.numeric(opts$threads))
registerDoSNOW(cl)

result <- foreach(i=1:length(query.names), .combine = bind_rows,
                  .packages = c("dplyr", "treeio", "tidytree", "phytools")) %dopar% {
  
  ind <- query.names[i]
  setTxtProgressBar(pb,i)
  ## Read in data, filter to most likely placement(s) and make various useful formats of it
  query.place.data <- in.jplace@placements %>% filter(name==ind)  %>% filter(like_weight_ratio >= 0.5*max(like_weight_ratio))
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
  
  
  
  ## Load thresholds for levels

  if(length(maxDistPlacements)==1)## If only one placement, maxDistPlacements isn't named, so add name
    names(maxDistPlacements) <- query.place.data$node
  
  ## Find lowest justifiable taxonomic classification of placement 
  ## Using mode of distances across placements to be conservative, for now
  numLevels <- lapply(maxDistPlacements, FUN=function(x) names(which(cutoffs>x))) %>% lapply(length) %>% getmode %>% unlist
  if(numLevels != 0){
    assignmentLevels <- names(cutoffs[1:numLevels])
  }else
    assignmentLevels <- NULL



  assignment <- tree.w.placements.tib[query.place.data$node, assignmentLevels]

  
  
  #if there is a multifurcation then assign the taxonomy as its parent node
  if(length(unique(tree.w.placements.tib$parent[query.place.data$node]))==1 & length(query.place.data$node)>1)
    assignment <- tree.w.placements.tib[unique(tree.w.placements.tib$parent[query.place.data$node]), rev(assignmentLevels)]

  
  assignment$query.name <- ind
  
  ### this is the part of the code that becomes problematic b/c it returns nothing if the index node doesn;t have an assignment
  ### how to deal with this...
  ### I guess the best way to deal with this is to figure out what 
  
  #rm_dist_index <- which(is.na(assignment[,names(cutoffs[numLevels])]))
  #so assignment will have a number of rows equal to the number of index nodes 
  #the assignments fron that index node will only go to the mode of the maxdistances 
  #we then want to return the maxdist so that it represents what is actually be assigned 
  #maybe we just change the name of the column from maxDist to median_distance
  #which represents the median distance of all potentially placements from the furthest child tip
  #this might be more interpretable 
  #in cases where there are multiple it should often represent the mode better
  #the other option is to select the nodes that are the mode cutoff and then return that
  #optionally we could keep the maxdist and retun the inidivaul distances for each assignment
  if(nrow(assignment)==1){
    assignment$maxDist <- max(maxDistPlacements)
  }else{
    maxDist_possible <- which(lapply(maxDistPlacements, FUN=function(x) names(which(cutoffs>x))) %>% lapply(length)==numLevels)
    assignment$maxDist <- max(maxDistPlacements[maxDist_possible])
  }

  
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

## potentially you could end up with cases where result is missing a column if no sequences are assigned at a particular level
taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
missing_levels <- which(!taxa_levels %in% colnames(result))

for(i in missing_levels){
  result[,taxa_levels[i]] <- NA
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

levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# for (i in 1:nrow(tax_parathaa)) {
#   if(grepl(";", tax_parathaa[i,"Kingdom"])){
#     tax_parathaa[i,levels[-1]] <- NA
#   }else if(grepl(";", tax_parathaa[i, "Phylum"])){
#     tax_parathaa[i, levels[-c(1,2)]] <- NA
#   }else if(grepl(";", tax_parathaa[i, "Class"])){
#     tax_parathaa[i, levels[-c(1,2,3)]] <- NA
#   }else if(grepl(";", tax_parathaa[i, "Order"])){
#     tax_parathaa[i, levels[-c(1,2,3,4)]] <- NA
#   }else if(grepl(";", tax_parathaa[i, "Family"])){
#     tax_parathaa[i, levels[-c(1,2,3,4,5)]] <- NA
#   }else if(grepl(";", tax_parathaa[i, "Genus"])){
#     tax_parathaa[i, levels[-c(1,2,3,4,5,6)]] <- NA
#   }
# }


if(delta>0){
  
  
  #make tree labels unique
  in.tree$label <- make.unique(in.tree$label)
  
  queries.w.species <- tax_parathaa %>% 
    filter(!is.na(Species)) %>% 
    pull(query.name)
  
  
  pb <- txtProgressBar(max=length(queries.w.species), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  pb_opts <- list(progress=progress)
  
  dists <- foreach(i=1:length(queries.w.species), .combine=bind_rows,
                   .packages = c("stringr", "castor")) %dopar% {
    
    #setTxtProgressBar(pb,i)
    query <- queries.w.species[i]
    temp <- nearest.neighbor.distances(tax.df=tax_parathaa, 
                                       placement.object=in.jplace, 
                                       reference.tree=in.tree, 
                                       max.radius=0.2,
                                       query = query)
    return(temp)
  }

##write.table(dists, file=file.path(opts$o, "distances.tsv"), sep='\t', quote=F, row.names=F)
  tax_parathaa <- nearest.neighbor.revisions(tax.df=tax_parathaa, 
                                             distances=dists, 
                                             radius=delta)
}
tax_parathaa$maxDist <- result$maxDist[match(tax_parathaa$query.name, result$query.name)]

write.table(tax_parathaa, file=file.path(opts$o, "taxonomic_assignments.tsv"),sep = '\t', quote = F, row.names=F)


