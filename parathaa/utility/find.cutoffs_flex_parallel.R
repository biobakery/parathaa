#!/usr/bin/env Rscript

### This script finds optimal cutoffs for each taxonomic level for a given tree, and plots their error scores

require(docopt)
'Usage:
   find.cutoffs.R [-d <naming file> -o <output> -n <tree> --wt1 <sweight> --wt2 <mweight> --threads <threads> --util1 <spec.editor_PATH> --util2 <cal.error_PATH>]

Options:
   -d naming file [default: input/taxmap_slv_ssu_ref_138.1.txt]
   -o model output directory [default: output/20230707_weight3Err10]
   -n tree [default: output/20230707_weight3Err10/region_specific.tree]
   --wt1 over-split penalty weight [default: 1]
   --wt2 over-merge penalty weight [default: 1]
   --threads number of threads to run in parallel [default: 1]
   --util1 path to utility file 1 [default: utility/SILVA.species.editor.dev.R]
   --util2 path to utility file 2 [default: utility/calc.error.scores.R]

 ]' -> doc

opts <- docopt(doc)


library(reshape2)
library(ggtree)
library(treeio)
library(tidyr)
library(dplyr)
library(ape)
library(ggplot2)
library(TDbook)
library(castor)
source(opts$util1)
source(opts$util2)
library(doSNOW)

#prevents plots from auto saving when called from Rscript
if(!interactive()) pdf(NULL)

## Bring in taxonomy file
inFileTaxdata <- opts$d

## Tree made from database trimmed to region
in.tree <- read.newick(opts$n)
in.tree.data <- as_tibble(in.tree)

## The above gives a warning because only the tip nodes have primary accessions and arbIDs.
## I suppressed it so that the warning message doesn't reach the user.


## Taxonomy of database from SILVA (primary accession, start, stop)

# I wonder if using a data.table/tibbly/big.data (pckge) would increase performance here?
# (below three lines run fairly slowly...)

taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="", check.names=FALSE)

isSILVA=FALSE

# the IDs that come with the file specific ID to the sequence but they also have start and stop for the sequence!
if("start" %in% colnames(taxdata)){
  isSILVA <- TRUE
  taxdata <- taxdata %>%
    unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
  
  suppressWarnings({
    taxdata <- taxdata %>%
      separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
      dplyr::rename(Species = organism_name)
  })
}else{
  taxdata$path <- gsub("\\|t__.*", "", taxdata$path)
  taxdata$path <- gsub("[a-z]__", "", taxdata$path)
  taxdata <- taxdata %>% separate(col=path, into=c("Kingdom", "Phylum", "Class", 
                                                   "Order", "Family", "Genus", "Species"), sep="\\|")
  
}

if(grepl("\\|M", in.tree.data$label[1])){
  suppressWarnings({
    in.tree.data <- in.tree.data %>%
      separate(col=label, into=c("arbID", "primaryAccession"), remove=F, sep="\\|")
  })
}else if(isSILVA){
  suppressWarnings({
    in.tree.data <- in.tree.data %>%
      separate(col=label, into=c("primaryAccession", "arbID"), remove=F, sep="\\.")
  })
}else{
  in.tree.data$primaryAccession <- in.tree.data$label
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
  in.tree.data <- SILVA.species.editor(in.tree.data)
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
cutoffs<-c(seq(0.001, 0.009, by=0.001), seq(0.01, 0.5, by=0.01), seq(0.55, 1.5, 0.05))


#Initialize
outputScores <- list()
outputNs <- list()
outputCounts <- list()
inputData <- in.tree.data
inputData$maxDists <- NA

# load in parallel back end
cl <- makeCluster(as.numeric(opts$threads))
registerDoSNOW(cl)
iterations <- length(which(inputData$isTip==F))
pb <- txtProgressBar(max=iterations, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
pb_opts <- list(progress=progress)
#eventually we should write code that doesn't run it in parallel with "1" worker if the threads is set to 1...
#get the maxdist of each internal node and the its tips

internal_node_stats <- foreach(i=1:length(which(inputData$isTip==F))) %dopar% {
  tmp_list <- list()
  
  intNode <- which(inputData$isTip==F)[i]
 #print(paste("Node:", intNode))
  
  tmp_node <- intNode - length(which(inputData$isTip==T))
  sub_tree <- castor::get_subtree_at_node(treeio::as.phylo(inputData), tmp_node)
  maxDist <- castor::find_farthest_tip_pair(sub_tree$subtree)$distance
  
  tmp_tips <- sub_tree$new2old_tip
  ch <- inputData[which(inputData$node %in% tmp_tips),]
  
  tmp_list[["maxDist"]] <- maxDist
  tmp_list[["tips"]] <- ch
  
  return(tmp_list)
}
save(internal_node_stats, file= file.path(opts$o, "internal_node_stats.RData")) ## used in next step
#set internal node names back on each cutoff distance
iterations <- length(cutoffs)
pb <- txtProgressBar(max=iterations, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
pb_opts <- list(progress=progress)


resultData <- foreach(i=1:length(cutoffs)) %dopar% {
  cut1 <- cutoffs[i]
  tmp_rez <- list()
  tmp_rez[[paste0("tax",cut1)]] <- inputData
  
  for(j in 1:length(internal_node_stats)){
    
    
    #map back to original node
    intNode <- which(inputData$isTip==F)[j]
    
    maxDist <- internal_node_stats[[j]][[1]]
    ch <- internal_node_stats[[j]][[2]]
    
    if(maxDist < cut1){
      if(length(table(ch[["Kingdom"]]))!=0)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Kingdom"] <- names(table(ch[["Kingdom"]]))[which(table(ch[["Kingdom"]])==max(table(ch[["Kingdom"]])))][[1]]
      if(length(table(ch[["Phylum"]]))!=0)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Phylum"] <- names(table(ch[["Phylum"]]))[which(table(ch[["Phylum"]])==max(table(ch[["Phylum"]])))][[1]]
      if(length(table(ch[["Class"]]))!=0)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Class"] <- names(table(ch[["Class"]]))[which(table(ch[["Class"]])==max(table(ch[["Class"]])))][[1]]
      if(length(table(ch[["Order"]]))!=0)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Order"] <- names(table(ch[["Order"]]))[which(table(ch[["Order"]])==max(table(ch[["Order"]])))][[1]]
      if(length(table(ch[["Family"]]))!=0)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Family"] <- names(table(ch[["Family"]]))[which(table(ch[["Family"]])==max(table(ch[["Family"]])))][[1]]
      if(length(table(ch[["Genus"]]))!=0)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Genus"] <- names(table(ch[["Genus"]]))[which(table(ch[["Genus"]])==max(table(ch[["Genus"]])))][[1]]
      if(length(table(ch[["Species"]]))!=0)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Species"] <- names(table(ch[["Species"]]))[which(table(ch[["Species"]])==max(table(ch[["Species"]])))][[1]]
    }
  }
  return(tmp_rez[[paste0("tax",cut1)]])
}
names(resultData) <- paste0("tax",cutoffs)
save(resultData, file= file.path(opts$o, "resultData.RData")) ## used in next step
# get the scores for each cutoff

iterations <- length(cutoffs)
pb <- txtProgressBar(max=iterations, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
pb_opts <- list(progress=progress)


outputScores <- foreach(i=1:length(cutoffs), .packages = c("dplyr", "treeio", "tidytree", "castor"), 
                        .combine = c) %dopar% {
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
  for(level in c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")){
    errs <- calc.error.scores(tempData, level, wt1=as.numeric(opts$wt1), wt2=as.numeric(opts$wt2))
    if(is.null(errs)){
      outputScores[[level]][as.character(cut1)] <- NA
    }else{
      outputScores[[level]][as.character(cut1)] <- errs[["scores"]] 
    }
  }
  ret_list <- list(outputScores)
  return(ret_list)
}

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
plotData2$Level <- factor(plotData2$Level, levels=c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", 
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

