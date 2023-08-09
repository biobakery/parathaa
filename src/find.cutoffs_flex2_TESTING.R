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
library(reshape2)

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
library(doParallel)


## Bring in taxonomy file
inFileTaxdata <- opts$d

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
cl <- makeCluster(as.numeric(opts$threads), outfile="~/Desktop/workers_test.txt")
registerDoParallel(cl)


#This loop needs to be re-written if we want it to run in parallel...

#something to do here...
start_time1 <- Sys.time()
foreach(i=1:length(which(inputData$isTip==F))) %do% {
  #progress bar?
  intNode <- which(inputData$isTip==F)[i]
  print(paste("Node:", intNode))
  
  tmp_node <- intNode - length(which(inputData$isTip==T))
  sub_tree <- castor::get_subtree_at_node(treeio::as.phylo(inputData), tmp_node)
  maxDist <- castor::find_farthest_tip_pair(sub_tree$subtree)$distance

  tmp_tips <- sub_tree$new2old_tip
  ch <- inputData[which(inputData$node %in% tmp_tips),]
  
  for(cut1 in cutoffs){

   
    if(maxDist < cut1) {
      #check if the cutoff is less than the maximum distance
      #if it is then we check if there are nodes underneath it assigned to X level
      
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
      if(length(table(ch[["Species"]]))!=0 & cut1 < 0.05)
        resultData[[paste0("tax",cut1)]][intNode, "Species"] <- names(table(ch[["Species"]]))[which(table(ch[["Species"]])==max(table(ch[["Species"]])))][[1]]
    }
  }
}
end_time1 <- Sys.time()
### this loop takes forever to run...
### running this in parallel is certainly possible..
#saveRDS(resultData, "~/Repos/Parathaa2_OP3/BB3/results_from_maxdists.RDS")



## update this loop to run in parallel...
start_time2 <- Sys.time()
outputScores <- foreach(i=1:length(cutoffs), .packages = c("dplyr", "treeio", "tidytree", "castor"), .combine = c) %dopar% {
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
  outputCounts <- list()
  outputScores <- list()
  for(level in c("Phylum", "Class", "Order", "Family", "Genus", "Species")){
    errs <- calc.error.scores(tempData, level, wt1=as.numeric(opts$wt1), wt2=as.numeric(opts$wt2))
    #outputCounts[[level]][[as.character(cut1)]] <- errs[["counts"]]
    outputScores[[level]][as.character(cut1)] <- errs[["scores"]]
  }
  ret_list <- list(outputScores)
  return(ret_list)
}
end_time2 <- Sys.time()

tmp_df <- do.call(cbind, outputScores)
colnames(tmp_df) <- names(unlist(tmp_df[1,]))
rnames <- rownames(tmp_df)

tmp_df <- apply(tmp_df, 2, function(x) unlist(x))
rownames(tmp_df) <- rnames
tmp_df_melt <- reshape2::melt(tmp_df)
colnames(tmp_df_melt) <- c("Level", "Threshold", "Scores")
## Create Plot of Error Scores for each threshold at each level
nseqs <- nrow(inputData %>% filter(isTip==TRUE))
tmp_df_melt$Scores <- tmp_df_melt$Scores/nseqs

plotData2 <- tmp_df_melt

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

stopCluster(cl)

