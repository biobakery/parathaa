#!/usr/bin/env Rscript

### This script finds optimal cutoffs for each taxonomic level for a given tree, and plots their error scores
if(!interactive()) pdf(NULL)

require(docopt)
'Usage:
   find.cutoffs.R [-d <naming file> -o <output> -n <tree> --wt1 <sweight> --wt2 <mweight> --bError <error_rate> --bThreshold <threshold>]

Options:
   -d naming file [default: input/taxmap_slv_ssu_ref_138.1.txt]
   -o model output directory [default: output/20230707_weight3Err10]
   -n tree [default: output/20230707_weight3Err10/region_specific.tree]
   --wt1 over-split penalty weight [default: 1]
   --wt2 over-merge penalty weight [default: 1]
   --bError binomial error rate [default: 1]
   --bThreshold binomial threshold [default: 1]
   --threads number of threads to run in parallel [default: 1]

 ]' -> doc

#To add as arguments: binom error params

opts <- docopt(doc)


# If we wanted we could add a loop that tries to install these packages within R?
# Although it might be nice to have this as a seperate R script 
# that runs on anadama workflow start up
library(logging)

# R logging example 
loginfo("Performing analysis data", logger="")

library(ggtree)
library(treeio)
library(tidyr)
library(dplyr)
library(ape)
#library(ggimage)
library(ggplot2)
library(TDbook)
library(castor)
source("src/SILVA.species.editor.R")
source("src/calc.error.scores.R")

## Bring in taxonomy file
inFileTaxdata <- opts$d

# inFileTaxdata <- "~/Repos/Parathaa2_OP3/BB3/gca2taxa.07022019.tsv"

## Tree made from database trimmed to region
in.tree <- read.newick(opts$n)
in.tree.data <- as_tibble(in.tree)


suppressWarnings({
  in.tree.data <- in.tree.data %>%
    separate(col=label, into=c("primaryAccession", "arbID"), remove=F, sep="\\.")
})
## The above gives a warning because only the tip nodes have primary accessions and arbIDs.
## I suppressed it so that the warning message doesn't reach the user.


## Taxonomy of database from SILVA (primary accession, start, stop)

# I wonder if using a data.table/tibbly/big.data (pckge) would increase performance here?
# (below three lines run fairly slowly...)

taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")

isSILVA=FALSE

# the IDs that come with the file specific ID to the sequence but they also have start and stop for the sequence!
if("start" %in% colnames(taxdata)){
  isSILVA <- TRUE
  taxdata <- taxdata %>%
    unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
  
  suppressWarnings({
    taxdata <- taxdata %>%
      separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
      dplyr::rename(Species = organism_name) %>%
      filter(Kingdom=="Bacteria") 
  })
}else{
  taxdata$path <- gsub("\\|t__.*", "", taxdata$path)
  taxdata$path <- gsub("[a-z]__", "", taxdata$path)
  taxdata <- taxdata %>% separate(col=path, into=c("Kingdom", "Phylum", "Class", 
                                                   "Order", "Family", "Genus", "Species"), sep="\\|") %>%
    filter(Kingdom=="Bacteria")
  
}

#in.tree.data is from the FAstTree... and taxdata is the silva taxonomy database
in.tree.data <- left_join(in.tree.data, taxdata, by="primaryAccession", multiple="all") %>%
  distinct(node, .keep_all=T) 

#some functions won't work unless its a specific type of class of data
class(in.tree.data) <- c("tbl_tree", class(in.tree.data))

#Check if each node is a tip or not!
in.tree.data$isTip <- isTip(in.tree.data, in.tree.data$node)

## Fix species names

# there are a lot of species names that are inconsistent and needed to be cleaned up!
if(isSILVA){
  in.tree.data <- SILVA.species.editor(in.tree.data, Task="assign_Tax")
}

in.tree.data <- in.tree.data %>% mutate(Kingdom = na_if(Kingdom, ""),
                                        Phylum = na_if(Phylum, ""),
                                        Class = na_if(Class, ""),
                                        Order = na_if(Order, ""),
                                        Family = na_if(Family, ""),
                                        Genus = na_if(Genus, ""),
                                        Species = na_if(Species, ""))  


# parameters from the binomModel
falseNegRate <- opts$bError
acceptableProb <- opts$bThreshold


#Broad range of cutoffs:
#looks at a series of cutoffs for the best "threshold"
#might be better to look minima and choose that and search around it rather than looking at all these...
cutoffs<-c(seq(0.001, 0.009, by=0.001), seq(0.01, 0.5, by=0.01), seq(0.55, 0.9, 0.05))


#Initialize
outputScores <- list()
outputNs <- list()
outputCounts <- list()
resultData <- list()
inputData <- in.tree.data
inputData$maxDists <- NA


## input data is the FastTree data and silva taxonomy bound together see above code
## for each threshold uses it as cutoff for that taxonomic level
## would be nice to have that dynamically chosen for each taxonomic level rather than using the same one...
for(cut1 in cutoffs){
  ## Note from Meg: This step is not currently memory-efficient, it stores a lot of the same data many times
  resultData[[paste0("tax",cut1)]] <- inputData
}

# Also a general comment! For loops and very slow in R I wonder if we can make use of other
# If we cannot vectorize this loop (using apply) maybe we can think about allow multi core
# loop using the "foreach" package?


col_data <- castor::collapse_tree_at_resolution(as.phylo(inputData),
  resolution = 0, shorten = F)

col_data_df <- as_tibble(col_data$tree) 
col_data_df$isTip <- isTip(col_data_df, col_data_df$node)

col_tips <- which(col_data_df$isTip==TRUE)

all_dists <- castor::get_all_pairwise_distances(col_data$tree, only_clades = col_tips)

internal_col_nodes <- which(col_data_df$isTip==FALSE) - length(col_data$tree$tip.label)

sub_trees <- get_subtrees_at_nodes(col_data$tree, internal_col_nodes)

## this is still very slow...
max_distance <- list()
for(col_intNode in internal_col_nodes){
  max_distance[[col_intNode]] <- max(all_dists[sub_trees$new2old_tips[[col_intNode]], sub_trees$new2old_tips[[col_intNode]]])
}




## go through each node that is not a tip

## okay this part is now extremely slow...
## how would we go about speeding it up?

## adding in parallel computation?
cl <- makeCluster(opts$threads, outfile="~/Desktop/workers_test.txt")
registerDoParallel(cl)



foreach(i=1:length(which(inputData$isTip==F))) %dopar% {
  #progress bar?
  intNode <- which(inputData$isTip==F)[i]
  print(paste("Node:", intNode))

  
#grab tips
  ch <- tidytree::offspring(inputData, intNode, tiponly=T)
  
  # cutoffs represent distances on the tree from tip to tip
  
  # do a correlation plot from distances to ANI would be nice!
  intNode_map <- col_data$old2new_clade[intNode]
  intNode_map

  if(intNode_map <= length(col_data$tree$tip.label)){
    maxDist <- 0
  }else{
    maxDist <- max_distance[[intNode_map - length(col_data$tree$tip.label)]]
  }
  
  for(cut1 in cutoffs){

   
    if(maxDist < cut1) {
      #check if the cutoff is less than the maximum distance
      #if it is then we check if there are nodes underneath it assigned to X level
      message(cut1)
      
      if(length(table(ch[["Kingdom"]]))!=0)
        #This looks like a crtical part of the code
        # We are making a new list item for the child node
        # We grab names of all the taxa at that level 
        # We assign the maxiumum name
        resultData[[paste0("tax",cut1)]][intNode, "Kingdom"] <- names(table(ch[["Kingdom"]]))[which(table(ch[["Kingdom"]])==max(table(ch[["Kingdom"]])))][[1]]
      if(length(table(ch[["Phylum"]]))!=0)
        resultData[[paste0("tax",cut1)]][intNode, "Phylum"] <- names(table(ch[["Phylum"]]))[which(table(ch[["Phylum"]])==max(table(ch[["Phylum"]])))][[1]]
      if(length(table(ch[["Class"]]))!=0)
        resultData[[paste0("tax",cut1)]][intNode, "Class"] <- names(table(ch[["Class"]]))[which(table(ch[["Class"]])==max(table(ch[["Class"]])))][[1]]
      if(length(table(ch[["Order"]]))!=0)
        resultData[[paste0("tax",cut1)]][intNode, "Order"] <- names(table(ch[["Order"]]))[which(table(ch[["Order"]])==max(table(ch[["Order"]])))][[1]]
      if(length(table(ch[["Family"]]))!=0)
        resultData[[paste0("tax",cut1)]][intNode, "Family"] <- names(table(ch[["Family"]]))[which(table(ch[["Family"]])==max(table(ch[["Family"]])))][[1]]
      if(length(table(ch[["Genus"]]))!=0)
        resultData[[paste0("tax",cut1)]][intNode, "Genus"] <- names(table(ch[["Genus"]]))[which(table(ch[["Genus"]])==max(table(ch[["Genus"]])))][[1]]
      if(length(table(ch[["Species"]]))!=0)
        resultData[[paste0("tax",cut1)]][intNode, "Species"] <- names(table(ch[["Species"]]))[which(table(ch[["Species"]])==max(table(ch[["Species"]])))][[1]]
    }
  }
}

### this loop takes forever to run...
### running this in parallel is certainly possible..


foreach(i=1:length(cutoffs)) %dopar% {
  cut1 <- cutoffs[[i]]
  print(cut1)
  tempData <- resultData[[paste0("tax",cut1)]]
  
### Define "Genus" as a clade defined by a Genus-named node with a non-Genus-named parent node, etc
  
  #if parent node is not kingdom assigned and this node is then its a kingdom node
  tempData$isKingdomNode <- is.na(tempData$Kingdom[tempData$parent]) & !is.na(tempData$Kingdom)
  tempData$isPhylumNode <- is.na(tempData$Phylum[tempData$parent]) & !is.na(tempData$Phylum)
  tempData$isClassNode <- is.na(tempData$Class[tempData$parent]) & !is.na(tempData$Class)
  tempData$isOrderNode <- is.na(tempData$Order[tempData$parent]) & !is.na(tempData$Order)
  tempData$isFamilyNode <- is.na(tempData$Family[tempData$parent]) & !is.na(tempData$Family)
  tempData$isGenusNode <- is.na(tempData$Genus[tempData$parent]) & !is.na(tempData$Genus)
  tempData$isSpeciesNode <- is.na(tempData$Species[tempData$parent]) & !is.na(tempData$Species)
  
  
  ## Calculate Error Scores for the given cutoff at each level
  
  for(level in c("Phylum", "Class", "Order", "Family", "Genus", "Species")){
    errs <- calc.error.scores(tempData, level, wt1=as.numeric(opts$wt1), wt2=as.numeric(opts$wt2))
    outputCounts[[level]][[as.character(cut1)]] <- errs[["counts"]]
    outputScores[[level]][as.character(cut1)] <- errs[["scores"]]
  }
  
}


## Create Plot of Error Scores for each threshold at each level
nseqs <- nrow(inputData %>% filter(isTip==TRUE))
plotData2 <- data.frame("Scores"=c(outputScores[["Kingdom"]]/nseqs,
                                   outputScores[["Phylum"]]/nseqs,
                                   outputScores[["Class"]]/nseqs,
                                   outputScores[["Order"]]/nseqs,
                                   outputScores[["Family"]]/nseqs,
                                   outputScores[["Genus"]]/nseqs, 
                                   outputScores[["Species"]]/nseqs),
                        "Level"=c(rep("Kingdom", length(outputScores[["Kingdom"]])), 
                                  rep("Phylum", length(outputScores[["Phylum"]])),
                                  rep("Class", length(outputScores[["Class"]])),
                                  rep("Order", length(outputScores[["Order"]])),
                                  rep("Family", length(outputScores[["Family"]])),
                                  rep("Genus", length(outputScores[["Genus"]])),
                                  rep("Species", length(outputScores[["Species"]]))),
                        "Threshold"=as.numeric(c(names(outputScores[["Kingdom"]]),
                                                 names(outputScores[["Phylum"]]),
                                                 names(outputScores[["Class"]]),
                                                 names(outputScores[["Order"]]),
                                                 names(outputScores[["Family"]]),
                                                 names(outputScores[["Genus"]]),
                                                 names(outputScores[["Species"]])))
)

# decide whatever cutoff gives the best low score!
plotData2$Level <- factor(plotData2$Level, levels=c("Phylum", "Class", "Order", "Family", "Genus", 
                                                    "Species"))
mins<- plotData2 %>% 
  group_by(Level) %>% 
  summarize(minScores = min(Scores))
plotData2 <- plotData2 %>% 
  left_join(mins) %>%
  group_by(Level) %>%
  mutate(minThreshold = Threshold[minScores==Scores][1])

ggplot(plotData2 , aes(x=Threshold, y=Scores, color=Level)) + geom_point() + 
  geom_vline(aes(xintercept = minThreshold, color=Level, linetype=Level)) + 
  theme(text = element_text(size = 14)) 

save(plotData2, file = file.path(opts$o, "optimal_scores.RData")) ## This is used in future steps
ggsave(filename = file.path(opts$o, "optimal_scores.png"), height = 4, width=5, units = "in")

stopCluster(cl = cl)
