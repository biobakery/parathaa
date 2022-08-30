#!/usr/bin/env Rscript
require(docopt)
'Usage:
   analysis.R [-d <naming file> -o <output> -b <blast result>]

Options:
   -d naming file [default: ~/Downloads/silvaCodes.txt]
   -o model output [default: results.csv]
   -b blast output [default: output/vsearch_results.txt]

 ]' -> doc

opts <- docopt(doc)

library(logging)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
library(overlapping)
library(gplots)
library(data.table)
source("src/plot.within.between.R")
source("src/heatmap.clade.R")
separateNames <- TRUE ## For Lea this should be FALSE, there is no separate naming file
taxonomySource <- "GTDB"

# R logging example 
loginfo("Performing analysis data", logger="")

library(ggtree)
library(treeio)
library(tidyr)
library(dplyr)
library(ape)
library(ggimage)
library(TDbook)

##Input and output files:

## GTDB IDs with 16S count > 1 and rep of species
#GTDBkeeps <- read.csv(inFileGTDBkeeps)
#GTDB4$GTDBcode <- substring(GTDB4$Query, 4)
#test <- GTDB4 %>% filter(GTDBcode %in% GTDBkeeps$organism_name)

LOCAL <- TRUE
if(LOCAL==TRUE){
inFileNR99 <- "/Users/mis696/proj/16s-region-checker/input/SSURefNR99_1200_slv_138_1_nds.ntree"
inFileTaxdata <- "/Users/mis696/proj/16s-region-checker/input/taxmap_slv_ssu_ref_138.1.txt"
inFileGTDB4 <- "/Users/mis696/proj/16s-region-checker/output/SILVA_NR99_GTDB_taxonomy.RData"
}
if(LOCAL==FALSE){
  outFileBac <- "/n/home05/mishort/16s-region-checker/output/Bactree.svg"
  inFileNR99 <- "/n/home05/mishort/16s-region-checker/input/SSURefNR99_1200_slv_138_1_nds_removespace.ntree"
  inFileTaxdata <- "/n/home05/mishort/16s-region-checker/input/taxmap_slv_ssu_ref_138.1.txt"
  inFileGTDB4 <- "/n/home05/mishort/16s-region-checker/input/SILVA_NR99_GTDB_taxonomy.RData"
}

## NR99 tree from online ARB file converted to Newick format (ARB IDs only)
##NR99 <- read.newick("/Users/mis696/proj/16s-region-checker/input/SSURefNR99_1200_slv_138_1_edit.tree")
##NR99 <- read.newick(inFileNR99)


## FastTree output from Seed database (combined IDs)
seed.tree1 <- read.newick("/Users/mis696/proj/Phylogeny_Taxonomy/SILVA_seed.tree")
## FastTree output from seed database trimmed to V4 region (combined IDs)
seed.treev4 <- read.newick("/Users/mis696/proj/Phylogeny_Taxonomy/SILVA_seed_V4.tree")

##seed.treev4 <- tree_subset(seed.treev4, node=5195)
v4treedata <- as_tibble(seed.treev4)
v4treedata <- v4treedata %>%
  separate(col=label, into=c("primaryAccession", "arbID"), remove=F)

## Taxonomy of NR99 database from SILVA (primary accession, start, stop)
taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
taxdata <- taxdata %>%
  separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  rename(Species = organism_name) %>%
  filter(Kingdom=="Bacteria") 

v4treedata <- left_join(v4treedata, taxdata, by="primaryAccession") %>%
  distinct(node, .keep_all=T) 

class(v4treedata) <- c("tbl_tree", class(v4treedata))
v4treedata$isTip <- isTip(v4treedata, v4treedata$node)
#v4treedata <- v4treedata %>%
#  mutate(Kingdom = replace(Kingdom, is.na(Kingdom) & isTip, ""),
#         Phylum = replace(Phylum, is.na(Phylum)& isTip, ""),
#         Class = replace(Class, is.na(Class)& isTip, ""),
#         Order = replace(Order, is.na(Order)& isTip, ""),
#         Family = replace(Family, is.na(Family)& isTip, ""),
#         Genus = replace(Genus, is.na(Genus)& isTip, ""),
#         Species = replace(Species, is.na(Species)& isTip, ""),
#  )
v4treedata <- v4treedata %>% mutate(Kingdom = na_if(Kingdom, ""),
                                Phylum = na_if(Phylum, ""),
                                Class = na_if(Class, ""),
                                Order = na_if(Order, ""),
                                Family = na_if(Family, ""),
                                Genus = na_if(Genus, ""),
                                Species = na_if(Species, ""))  


if(FALSE){
## Results from vsearch of GTDB 16S sequences against NR99 db--~13k hits (primaryAccession.start.stop)
##load(file="~/proj/16s-region-checker/output/SILVA_GTDB_taxonomy.R") # GTDB4
load(inFileGTDB4) # GTDB4
## Separating out IDs
GTDB4 <- GTDB4 %>%
  separate(col=Match, into=c("primaryAccession", "start", "stop"), remove=F)## %>%
  ##distinct(Match, .keep_all=T)

taxdata <-  left_join(taxdata, GTDB4, by=c("AccID"="Match")) 
taxdata <- taxdata %>% mutate("GTDBmapped"=!is.na(Kingdom_GTDB))
taxdata <- taxdata %>% select(!primaryAccession.y) %>%
  rename(primaryAccession=primaryAccession.x)
}


Task <- "assign_Tax" # "assign_Tax" "find_cutoffs"
binomErrModel <- TRUE
falseNegRate <- 0.05
acceptableProb <- 0.20


#Broad range of cutoffs:
if(Task=="find_cutoffs")
  cutoffs<-c(seq(0.01, 0.5, by=0.01), seq(0.55, 0.9, 0.05))
#"Best" cutoffs:
if(Task=="assign_Tax")
  cutoffs<-c("Species"=0.01, "Genus"=0.04, "Family"=0.11, "Order"=0.17,  "Class"=0.3, "Phylum"=0.37)

#Initialize
outputScores <- list()
outputNs <- list()
outputCounts <- list()

resultData <- list()
inputData <- v4treedata
#inputData <- test
inputData$maxDists <- NA

if(Task=="find_cutoffs"){
  for(cut1 in cutoffs){
    resultData[[paste0("tax",cut1)]] <- inputData
  }
}
if(Task=="assign_Tax")
  resultData[["tax_bestcuts"]] <- inputData
startTime <- Sys.time()
for (intNode in which(inputData$isTip==F)){
  print(paste("Node:", intNode))
  ch <- offspring(inputData, intNode, tiponly=T)
  tre <- tree_subset(as.treedata(inputData), intNode, levels_back = 0)
  maxDist <- max(cophenetic(as.phylo(tre)))
  
  if(Task=="find_cutoffs"){
    for(cut1 in cutoffs){
      if(maxDist < cut1) {
        if(length(table(ch[["Kingdom"]]))!=0)
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
  } else if(Task=="assign_Tax" & binomErrModel==TRUE){
    resultData[["tax_bestcuts"]][intNode, "maxDists"] <- maxDist
    if(length(table(ch[["Kingdom"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Kingdom"] <- names(table(ch[["Kingdom"]]))[which(table(ch[["Kingdom"]])==max(table(ch[["Kingdom"]])))][[1]]
    if(maxDist < cutoffs["Phylum"] & length(table(ch[["Phylum"]]))!=0){
      if(pbinom(sum(table(ch[["Phylum"]]))-max(table(ch[["Phylum"]]))-1,
                sum(table(ch[["Phylum"]])),falseNegRate, lower.tail = F) > acceptableProb){
        resultData[["tax_bestcuts"]][intNode, "Phylum"] <- names(table(ch[["Phylum"]]))[which(table(ch[["Phylum"]])==max(table(ch[["Phylum"]])))][[1]]
      } else {resultData[["tax_bestcuts"]][intNode, "Phylum"] <- paste(
        names(which(table(ch[["Phylum"]]) >
                      qbinom(acceptableProb , sum(table(ch[["Phylum"]])),falseNegRate, lower.tail = F))), collapse=";")}
    }
    if(maxDist < cutoffs["Class"] & length(table(ch[["Class"]]))!=0){
      if(pbinom(sum(table(ch[["Class"]]))-max(table(ch[["Class"]]))-1,
                sum(table(ch[["Class"]])),falseNegRate, lower.tail = F) > acceptableProb){
        resultData[["tax_bestcuts"]][intNode, "Class"] <- names(table(ch[["Class"]]))[which(table(ch[["Class"]])==max(table(ch[["Class"]])))][[1]]
      } else {resultData[["tax_bestcuts"]][intNode, "Class"] <- paste(
        names(which(table(ch[["Class"]]) >
                      qbinom(acceptableProb , sum(table(ch[["Class"]])),falseNegRate, lower.tail = F))), collapse=";")}
    }
    if(maxDist < cutoffs["Order"] & length(table(ch[["Order"]]))!=0){
      if(pbinom(sum(table(ch[["Order"]]))-max(table(ch[["Order"]]))-1,
                sum(table(ch[["Order"]])),falseNegRate, lower.tail = F) > acceptableProb){
        resultData[["tax_bestcuts"]][intNode, "Order"] <- names(table(ch[["Order"]]))[which(table(ch[["Order"]])==max(table(ch[["Order"]])))][[1]]
      } else {resultData[["tax_bestcuts"]][intNode, "Order"] <- paste(
        names(which(table(ch[["Order"]]) >
                      qbinom(acceptableProb , sum(table(ch[["Order"]])),falseNegRate, lower.tail = F))), collapse=";")}
    }
    if(maxDist < cutoffs["Family"] & length(table(ch[["Family"]]))!=0){
      if(pbinom(sum(table(ch[["Family"]]))-max(table(ch[["Family"]]))-1,
                sum(table(ch[["Family"]])),falseNegRate, lower.tail = F) > acceptableProb){
        resultData[["tax_bestcuts"]][intNode, "Family"] <- names(table(ch[["Family"]]))[which(table(ch[["Family"]])==max(table(ch[["Family"]])))][[1]]
      } else {resultData[["tax_bestcuts"]][intNode, "Family"] <- paste(
        names(which(table(ch[["Family"]]) >
                      qbinom(acceptableProb , sum(table(ch[["Family"]])),falseNegRate, lower.tail = F))), collapse=";")}
    }
    if(maxDist < cutoffs["Genus"] & length(table(ch[["Genus"]]))!=0){
      if(pbinom(sum(table(ch[["Genus"]]))-max(table(ch[["Genus"]]))-1,
                sum(table(ch[["Genus"]])),falseNegRate, lower.tail = F) > acceptableProb){
        resultData[["tax_bestcuts"]][intNode, "Genus"] <- names(table(ch[["Genus"]]))[which(table(ch[["Genus"]])==max(table(ch[["Genus"]])))][[1]]
      } else {resultData[["tax_bestcuts"]][intNode, "Genus"] <- paste(
        names(which(table(ch[["Genus"]]) >
                      qbinom(acceptableProb , sum(table(ch[["Genus"]])),falseNegRate, lower.tail = F))), collapse=";")}
    }
    if(maxDist < cutoffs["Species"] & length(table(ch[["Species"]]))!=0){
      if(pbinom(sum(table(ch[["Species"]]))-max(table(ch[["Species"]]))-1,
                sum(table(ch[["Species"]])),falseNegRate, lower.tail = F) > acceptableProb){
        resultData[["tax_bestcuts"]][intNode, "Species"] <- names(table(ch[["Species"]]))[which(table(ch[["Species"]])==max(table(ch[["Species"]])))][[1]]
      } else {resultData[["tax_bestcuts"]][intNode, "Species"] <- paste(
        names(which(table(ch[["Species"]]) >
                      qbinom(acceptableProb , sum(table(ch[["Species"]])),falseNegRate, lower.tail = F))), collapse=";")}
    }

  } else if(Task=="assign_Tax" & binomErrModel==FALSE){
    resultData[["tax_bestcuts"]][intNode, "maxDists"] <- maxDist
    if(length(table(ch[["Kingdom"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Kingdom"] <- names(table(ch[["Kingdom"]]))[which(table(ch[["Kingdom"]])==max(table(ch[["Kingdom"]])))][[1]]
    if(maxDist < cutoffs["Phylum"] & length(table(ch[["Phylum"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Phylum"] <- names(table(ch[["Phylum"]]))[which(table(ch[["Phylum"]])==max(table(ch[["Phylum"]])))][[1]]
    if(maxDist < cutoffs["Class"] & length(table(ch[["Class"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Class"] <- names(table(ch[["Class"]]))[which(table(ch[["Class"]])==max(table(ch[["Class"]])))][[1]]
    if(maxDist < cutoffs["Order"] & length(table(ch[["Order"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Order"] <- names(table(ch[["Order"]]))[which(table(ch[["Order"]])==max(table(ch[["Order"]])))][[1]]
    if(maxDist < cutoffs["Family"] &length(table(ch[["Family"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Family"] <- names(table(ch[["Family"]]))[which(table(ch[["Family"]])==max(table(ch[["Family"]])))][[1]]
    if(maxDist < cutoffs["Genus"] & length(table(ch[["Genus"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Genus"] <- names(table(ch[["Genus"]]))[which(table(ch[["Genus"]])==max(table(ch[["Genus"]])))][[1]]
    if(maxDist < cutoffs["Species"] &length(table(ch[["Species"]]))!=0)
      resultData[["tax_bestcuts"]][intNode, "Species"] <- names(table(ch[["Species"]]))[which(table(ch[["Species"]])==max(table(ch[["Species"]])))][[1]]
  }
    
}

if(Task=="assign_Tax" & binomErrModel==T){
  ### Define "Genus" as a clade defined by a Genus-named node with a non-Genus-named parent node, etc
  resultData[["tax_bestcuts"]]$isPhylumNode <- is.na(resultData[["tax_bestcuts"]]$Phylum[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Phylum)
  resultData[["tax_bestcuts"]]$isClassNode <- is.na(resultData[["tax_bestcuts"]]$Class[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Class)
  resultData[["tax_bestcuts"]]$isOrderNode <- is.na(resultData[["tax_bestcuts"]]$Order[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Order)
  resultData[["tax_bestcuts"]]$isFamilyNode <- is.na(resultData[["tax_bestcuts"]]$Family[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Family)
  resultData[["tax_bestcuts"]]$isGenusNode <- is.na(resultData[["tax_bestcuts"]]$Genus[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Genus)
  resultData[["tax_bestcuts"]]$isSpeciesNode <- is.na(resultData[["tax_bestcuts"]]$Species[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Species)
}

elapsedTime <- Sys.time() - startTime
print(paste("Time elapsed:", elapsedTime, "minutes"))

for(cut1 in cutoffs){
  print(cut1)
  if(Task=="find_cutoffs")
    tempData <- resultData[[paste0("tax",cut1)]]
  if(Task=="assign_Tax")
    tempData <- resultData[["tax_bestcuts"]]
  ### Define "Genus" as a clade defined by a Genus-named node with a non-Genus-named parent node, etc
  tempData$isKingdomNode <- is.na(tempData$Kingdom[tempData$parent]) & !is.na(tempData$Kingdom)
  tempData$isPhylumNode <- is.na(tempData$Phylum[tempData$parent]) & !is.na(tempData$Phylum)
  tempData$isClassNode <- is.na(tempData$Class[tempData$parent]) & !is.na(tempData$Class)
  tempData$isOrderNode <- is.na(tempData$Order[tempData$parent]) & !is.na(tempData$Order)
  tempData$isFamilyNode <- is.na(tempData$Family[tempData$parent]) & !is.na(tempData$Family)
  tempData$isGenusNode <- is.na(tempData$Genus[tempData$parent]) & !is.na(tempData$Genus)
  tempData$isSpeciesNode <- is.na(tempData$Species[tempData$parent]) & !is.na(tempData$Species)
  
  ## Calculate Phylum Scores
  pNodes <- tempData %>% filter(isPhylumNode==T)
  ## First score to minimize: number of Phylum represented multiple times
  ## Second score to minimize: number of non-Phylum members clustered within a Phylum
  phyla <- offspring(tempData, tempData$node[which(tempData$isPhylumNode==T)], tiponly=T, self_include=T)
  phylumTables <- lapply(phyla, FUN= function(x) table(x[["Phylum"]]))
  phylumNames <- tempData %>% filter(node %in% names(phylumTables)) %>% select(Phylum, node)
  phylumNames$score2count <-  unlist(lapply(phylumTables, FUN=function(x) sum(x)-max(x)))
  phylumscore2s <- phylumNames %>% group_by(Phylum) %>% summarize(score2count = sum(score2count))
  phylumscore2s$score1count <- as.integer(table(pNodes$Phylum)-1)
  pscore1 <- sum(phylumscore2s$score1count)
  pscore2 <- sum(phylumscore2s$score2count)
  pscoreSum <- pscore1 + pscore2
  phylumscore2s$scoreSum <- phylumscore2s$score1count + phylumscore2s$score2count
  if(Task=="find_cutoffs"){
    outputCounts[["Phylum"]][[as.character(cut1)]] <- phylumscore2s
    outputScores[["Phylum"]][as.character(cut1)] <- pscoreSum
    #outputNs[["Phylum"]][as.character(cut1)] <- (lapply(phyla,nrow))
  }
  if(Task=="assign_Tax" & binomErrModel==F){
    assignGT1cut <- quantile(phylumNames$score2count, 0.95)
    assignGT1 <- phylumNames$node[which(phylumNames$score2count>assignGT1cut)]
    tempData[assignGT1, "Phylum"] <- lapply(
      lapply(phyla[which(names(phyla) %in% assignGT1)], FUN=function(x) unique(x$Phylum)), 
      FUN=function(x) paste(x, collapse=';')) %>% unlist()
  }
  
  ## Calculate Class Scores
  cNodes <- tempData %>% filter(isClassNode==T)
  ## First score to minimize: number of Class represented multiple times
  ## Second score to minimize: number of non-Class members clustered within a Class
  classes <- offspring(tempData, tempData$node[which(tempData$isClassNode==T)], tiponly=T, self_include=T)
  classTables <- lapply(classes, FUN= function(x) table(x[["Class"]]))
  classNames <- tempData %>% filter(node %in% names(classTables)) %>% select(Class, node)
  classNames$score2count <-  unlist(lapply(classTables, FUN=function(x) sum(x)-max(x)))
  classscore2s <- classNames %>% group_by(Class) %>% summarize(score2count = sum(score2count))
  classscore2s$score1count <- as.integer(table(cNodes$Class)-1)
  cscore1 <- sum(classscore2s$score1count)
  cscore2 <- sum(classscore2s$score2count)
  cscoreSum <- cscore1 + cscore2
  classscore2s$scoreSum <- classscore2s$score1count + classscore2s$score2count
  if(Task=="find_cutoffs"){
    outputCounts[["Class"]][[as.character(cut1)]] <- classscore2s
    outputScores[["Class"]][as.character(cut1)] <- cscoreSum
    #outputNs[["Class"]][as.character(cut1)] <- (lapply(classes,nrow))
  }
  if(Task=="assign_Tax" & binomErrModel==F){
    assignGT1cut <- quantile(classNames$score2count, 0.95)
    assignGT1 <- classNames$node[which(classNames$score2count>assignGT1cut)]
    tempData[assignGT1, "Class"] <- lapply(
      lapply(classes[which(names(classes) %in% assignGT1)], FUN=function(x) unique(x$Class)), 
      FUN=function(x) paste(x, collapse=';')) %>% unlist()
  }
  
  ## Calculate Order Scores
  oNodes <- tempData %>% filter(isOrderNode==T)
  ## First score to minimize: number of Orders represented multiple times
  ## Second score to minimize: number of non-Order members clustered within a Order
  orders <- offspring(tempData, tempData$node[which(tempData$isOrderNode==T)], tiponly=T, self_include=T)
  orderTables <- lapply(orders, FUN= function(x) table(x[["Order"]]))
  orderNames <- tempData %>% filter(node %in% names(orderTables)) %>% select(Order, node)
  orderNames$score2count <-  unlist(lapply(orderTables, FUN=function(x) sum(x)-max(x)))
  orderscore2s <- orderNames %>% group_by(Order) %>% summarize(score2count = sum(score2count))
  orderscore2s$score1count <- as.integer(table(oNodes$Order)-1)
  oscore1 <- sum(orderscore2s$score1count) 
  oscore2 <- sum(orderscore2s$score2count)
  oscoreSum <- oscore1 + oscore2
  orderscore2s$scoreSum <- orderscore2s$score1count + orderscore2s$score2count
  if(Task=="find_cutoffs"){
    outputCounts[["Order"]][[as.character(cut1)]] <- orderscore2s
    outputScores[["Order"]][as.character(cut1)] <- oscoreSum
    #outputNs[["Order"]][as.character(cut1)] <- (lapply(orders,nrow))
  }
  if(Task=="assign_Tax" & binomErrModel==F){
    assignGT1cut <- quantile(orderNames$score2count, 0.95)
    assignGT1 <- orderNames$node[which(orderNames$score2count>assignGT1cut)]
    tempData[assignGT1, "Order"] <- lapply(
      lapply(orders[which(names(orders) %in% assignGT1)], FUN=function(x) unique(x$Order)), 
      FUN=function(x) paste(x, collapse=';')) %>% unlist()
  }
  
  ## Calculate Family Scores
  fNodes <- tempData %>% filter(isFamilyNode==T)
  ## First score to minimize: number of Family represented multiple times
  ## Second score to minimize: number of non-Family members clustered within a Family
  families <- offspring(tempData, tempData$node[which(tempData$isFamilyNode==T)], tiponly=T, self_include=T)
  familyTables <- lapply(families, FUN= function(x) table(x[["Family"]]))
  familyNames <- tempData %>% filter(node %in% names(familyTables)) %>% select(Family, node)
  familyNames$score2count <- unlist(lapply(familyTables, FUN=function(x) sum(x)-max(x)))
  familyscore2s <- familyNames %>% group_by(Family) %>% summarize(score2count = sum(score2count))
  familyscore2s$score1count <- as.integer(table(fNodes$Family)-1)
  fscore1 <- sum(familyscore2s$score1count)
  fscore2 <- sum(familyscore2s$score2count)
  fscoreSum <- fscore1 + fscore2
  familyscore2s$scoreSum <- familyscore2s$score1count + familyscore2s$score2count
  if(Task=="find_cutoffs"){
    outputCounts[["Family"]][[as.character(cut1)]] <- familyscore2s
    outputScores[["Family"]][as.character(cut1)] <- fscoreSum
    #outputNs[["Family"]][as.character(cut1)] <- (lapply(families,nrow))
  }
  if(Task=="assign_Tax" & binomErrModel==F){
    assignGT1cut <- quantile(familyNames$score2count, 0.95)
    assignGT1 <- familyNames$node[which(familyNames$score2count>assignGT1cut)]
    tempData[assignGT1, "Family"] <- lapply(
      lapply(families[which(names(families) %in% assignGT1)], FUN=function(x) unique(x$Family)), 
      FUN=function(x) paste(x, collapse=';')) %>% unlist()
  }
  
  ## Calculate Genus Scores
  gNodes <- tempData %>% filter(isGenusNode==T)
  ## First score to minimize: number of Genera represented multiple times
  ## Second score to minimize: number of non-Genus members clustered within a Genus
  genera <- offspring(tempData, tempData$node[which(tempData$isGenusNode==T)], tiponly=T, self_include=T)
  generaTables <- lapply(genera, FUN= function(x) table(x[["Genus"]]))
  generaNames <- tempData %>% filter(node %in% names(generaTables)) %>% select(Genus, node)
  generaNames$score2count <-  unlist(lapply(generaTables, FUN=function(x) sum(x)-max(x)))
  generascore2s <- generaNames %>% group_by(Genus) %>% summarize(score2count = sum(score2count))
  generascore2s$score1count <- as.integer(table(gNodes$Genus)-1)
  gscore1 <- sum(generascore2s$score1count)
  gscore2 <- sum(generascore2s$score2count)
  gscoreSum <- gscore1 + gscore2
  generascore2s$scoreSum <- generascore2s$score1count + generascore2s$score2count
  if(Task=="find_cutoffs"){
    outputCounts[["Genus"]][[as.character(cut1)]] <- generascore2s
    outputScores[["Genus"]][as.character(cut1)] <- gscoreSum
    #outputNs[["Genus"]][as.character(cut1)] <- (lapply(genera,nrow))
  }
  if(Task=="assign_Tax" & binomErrModel==F){
    assignGT1cut <- quantile(generaNames$score2count, 0.95)
    assignGT1 <- generaNames$node[which(generaNames$score2count>assignGT1cut)]
    tempData[assignGT1, "Genus"] <- lapply(
      lapply(genera[which(names(genera) %in% assignGT1)], FUN=function(x) unique(x$Genus)), 
      FUN=function(x) paste(x, collapse=';')) %>% unlist()
  }
  
  ## Calculate Species Scores
  sNodes <- tempData %>% filter(isSpeciesNode==T)
  ## First score to minimize: number of Species represented multiple times
  ## Second score to minimize: number of non-Species members clustered within a Species
  species <- offspring(tempData, tempData$node[which(tempData$isSpeciesNode==T)], tiponly=T, self_include=T)
  speciesTables <- lapply(species, FUN= function(x) table(x[["Species"]]))
  speciesNames <- tempData %>% filter(node %in% names(speciesTables)) %>% select(Species, node)
  speciesNames$score2count <-  unlist(lapply(speciesTables, FUN=function(x) sum(x)-max(x)))
  speciesNames$nSeqs <- unlist(lapply(species, nrow))
  speciesscore2s <- speciesNames %>% group_by(Species) %>% summarize(score2count = sum(score2count))
  speciesscore2s$score1count <- as.integer(table(sNodes$Species)-1)
  sscore1 <- sum(speciesscore2s$score1count)
  sscore2 <- sum(speciesscore2s$score2count)
  sscoreSum <- sscore1 + sscore2
  speciesscore2s$scoreSum <- speciesscore2s$score1count + speciesscore2s$score2count
  if(Task=="find_cutoffs"){
    outputCounts[["Species"]][[as.character(cut1)]] <- speciesscore2s
    outputScores[["Species"]][as.character(cut1)] <- sscoreSum
    #outputNs[["Species"]][as.character(cut1)] <- (lapply(species,nrow))
  }
  if(Task=="assign_Tax" & binomErrModel==F){
    assignGT1cut <- quantile(speciesNames$score2count, 0.95)
    assignGT1 <- speciesNames$node[which(speciesNames$score2count>assignGT1cut)]
    tempData[assignGT1, "Species"] <- lapply(
      lapply(species[which(names(species) %in% assignGT1)], FUN=function(x) unique(x$Species)), 
      FUN=function(x) paste(x, collapse=';')) %>% unlist()
    resultData[["tax_bestcuts"]] <- tempData
    break()
  }
}

Sys.time()

#Example tree
test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 7104, levels_back = 7)
ggtree(test1) + geom_tippoint(aes(color=Phylum, size=1.5)) + geom_nodepoint(aes(color=Phylum, shape=isPhylumNode, size=1.5))
ggsave(filename = "~/proj/16s-region-checker/output/sampleAmbigTree_20220510.png" , height=10, width=8, units="in")


test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 7269, levels_back = 2)
ggtree(test1) + geom_tippoint(aes(color=Family, size=1.5)) + geom_nodepoint(aes(color=Family, shape=isFamilyNode, size=1.5))+ geom_text(aes(label=label), hjust=-.1)
ggtree(test1) + geom_tippoint(aes(color=Genus, size=1.5)) + geom_nodepoint(aes(color=Genus, shape=isGenusNode, size=1.5))


test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 7104, levels_back = 1)
ggtree(test1) + geom_tippoint(aes(color=Phylum, size=1.5)) + geom_nodepoint(aes(color=Phylum, shape=isPhylumNode, size=1.5))

# Bacteroides
test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 6119, levels_back =15)
ggtree(test1) + geom_tippoint(aes(color=Phylum, size=1.5)) + geom_nodepoint(aes(color=Phylum, shape=isPhylumNode, size=1.5))
test1 <- tree_subset(as.treedata(resultDataBin$tax_bestcuts), 6119, levels_back = 15)
ggtree(test1) + geom_tippoint(aes(color=Phylum, size=1.5)) + geom_nodepoint(aes(color=Phylum, shape=isPhylumNode, size=1.5))

test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 8063, levels_back = 1)
ggtree(test1) + geom_tippoint(aes(color=Family, size=1.5)) + geom_nodepoint(aes(color=Family, shape=isFamilyNode, size=1.5))

test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 4764, levels_back = 10)
ggtree(test1) + geom_tippoint(aes(color=Order, size=1.5)) + geom_nodepoint(aes(color=Order, shape=isOrderNode, size=1.5)) 
ggtree(test1) + geom_tippoint(aes(color=Family, size=1.5)) + geom_nodepoint(aes(color=Family, shape=isFamilyNode, size=1.5))

ggtree(test1) + geom_tippoint(aes(color=Genus, size=1.5)) + geom_nodepoint(aes(color=Genus, shape=isGenusNode, size=1.5))

test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 4764, levels_back = 3)
ggtree(test1) + geom_tippoint(aes(color=Species, size=1.5)) + geom_nodepoint(aes(color=Species, shape=isSpeciesNode, size=1.5)) 

test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 5811, levels_back=0)
ggtree(test1)+ geom_tippoint(aes(color=Phylum, size=1.5)) + geom_nodepoint(aes(color=Phylum, shape=isPhylumNode, size=1.5))


test1 <- tree_subset(as.treedata(resultData$tax_bestcuts), 6119, levels_back =1)
ggtree(test1)+ geom_tippoint(aes(color=Species, size=1.5)) + geom_nodepoint(aes(color=Species, shape=isSpeciesNode, size=1.5)) + geom_text(aes(label=node), hjust=-.5)
ggtree(test1)+ geom_tippoint(aes(color=Genus, size=1.5)) + geom_nodepoint(aes(color=Genus, shape=isGenusNode, size=1.5)) + geom_text(aes(label=node), hjust=-.5)
ggtree(test1)+ geom_tippoint(aes(color=Family, size=1.5)) + geom_nodepoint(aes(color=Family, shape=isFamilyNode, size=1.5)) + geom_text(aes(label=node), hjust=-.5)
ggtree(test1)+ geom_tippoint(aes(color=Order, size=1.5)) + geom_nodepoint(aes(color=Order, shape=isOrderNode, size=1.5)) + geom_text(aes(label=node), hjust=-.5)
test2 <- tree_subset(seed.treev4, 6119, levels_back =1)
test3 <- multi2di(test2)
ape::write.tree(test3, file='~/proj/Phylogeny_Taxonomy/bacTree_small.nwk')
as_tibble(test3) %>% View()
write.table(as_tibble(test3)$label[isTip(test3, as_tibble(test3)$node)], quote = F, file = '~/proj/Phylogeny_Taxonomy/bacTreeNames_small.txt', row.names = F, col.names = F)
p = ggtree(test1)+ geom_tippoint(aes(color=Genus, size=1.5)) + geom_nodepoint(aes(color=Genus, shape=isGenusNode, size=1.5)) + geom_tiplab()


##Run for optimal Family cutoff and plot across families:
##Error score per family:
if(FALSE){
  pct95Fam <- quantile(fScores$scoreCount/nseqs, 0.95)
  png(filename="proj/16s-region-checker/output/FamilyErrors_fullLength_20211214.png")
  hist(familyScores$scoreCount/nseqs, xlab="Error Score for Family", main="")
  abline(v=pct95Fam, col="blue")
  text(x=0.004, y=150, paste("95th %ile:", round(pct95Fam,digits = 4 )), col="blue")
  dev.off()
}
if(FALSE){
  pct95Gen <- quantile(generaScores$scoreCount/nseqs, 0.95)
  png(filename="proj/16s-region-checker/output/GenusErrors_fullLength_20211214.png")
  hist(generaScores$scoreCount/nseqs, xlab="Error Score for Genera", main="")
  abline(v=pct95Gen, col="blue")
  text(x=0.005, y=150, paste("95th %ile:", round(pct95Gen,digits = 4 )), col="blue")
  dev.off()
}

if(FALSE){ #temp
## Optimal Thresholds Plot
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
  plotData2$Level <- factor(plotData2$Level, levels=c("Phylum", "Class", "Order", "Family", "Genus", 
                                                       "Species"))
mins<- plotData2 %>% 
  group_by(Level) %>% 
  summarize(minScores = min(Scores))
plotData2 <- plotData2 %>% 
  left_join(mins) %>%
  group_by(Level) %>%
  mutate(minThreshold = Threshold[minScores==Scores])
  
ggplot(plotData2 , aes(x=Threshold, y=Scores, color=Level)) + geom_point() + 
  geom_vline(aes(xintercept = minThreshold, color=Level, linetype=Level))

save(plotData2, file = "~/proj/16s-region-checker/output/scores_thresholdOnly20220315.R")
ggsave(filename="~/proj/16s-region-checker/output/scores_thresholdOnly_20220315.png")

## Optimal Thresholds Plot--full-length 16S gene
plotData3 <- data.frame("Scores"=c(outputScoresGeneraFull/nseqs,outputScoresFamilyFull/nseqs), 
                        "Level"=c(rep("Genus", length(outputScoresGeneraFull)), rep("Family", length(outputScoresFamilyFull))),
                        "Threshold"=as.numeric(c(names(outputScoresGeneraFull), names(outputScoresFamilyFull)))
)
ggplot(plotData3, aes(x=Threshold, y=Scores, color=Level)) + geom_point()
#ggsave(filename="~/proj/16s-region-checker/output/scores_thresholdOnly_fullLength20211214.png")
#save(plotData3, file = "~/proj/16s-region-checker/output/scores_thresholdOnly_fullLength20211214.R")


ggtree(as.treedata(v4data)) + geom_tippoint(aes(color=Genus, shape=isGenusNode)) + geom_nodepoint(aes(color=Genus, shape=isGenusNode)) 

#p.v4sub <- ggtree(seed.treev4.subset) %<+% taxdata + geom_tippoint(aes(color=Order))  + geom_tiplab() + geom_nodelab()
p.v4sub <-  ggtree(as.treedata(v4data))+ geom_tippoint(aes(color=Genus)) + geom_nodepoint(aes(color=Genus)) 
p.v4sub + layout_dendrogram()

max(cophenetic(as.phylo(v4subset.wdata)))

cophenetic(as.phylo(v4subset.wdata), check=T)

v4Int <- v4data %>% filter(isTip==F) 
table(is.na(v4Int$Genus))


node292 <- tree_subset(as.treedata(v4data), node="0.988", levels_back=0)
max(cophenetic(as.phylo(node292)))
ggtree(node292)+ geom_tippoint(aes(color=Family)) 

node292 <- tree_subset(as.treedata(v4data), node="0.292", levels_back=0)
max(cophenetic(as.phylo(node292)))
ggtree(node292)+ geom_tippoint(aes(color=Family)) 

node560 <- tree_subset(as.treedata(v4data), node="0.560", levels_back=0)
max(cophenetic(as.phylo(node560)))
ggtree(node560)+ geom_tippoint(aes(color=Family)) 
 

Bac <- fNodes %>% filter(Family=="Bacillaceae")


BacNodes <- tree_subset(as.treedata(v4data), node=1664, levels_back=8)
ggtree(BacNodes)+ geom_tippoint(aes(color=Family, shape=isFamilyNode)) +geom_nodepoint(aes(color=Family, shape=isFamilyNode))

## Tree of full-length sequences
p <- 
  ggtree(seed.tree1, branch.length='none', layout='circular') %<+% taxdata[,1:12] + geom_tippoint(aes(color=Phylum)) 
## Tree of v4 sequences
p.v4 <- 
  ggtree(seed.treev4, branch.length='none', layout='circular') %<+% taxdata + geom_tippoint(aes(color=Order))
p.v4
}#temp



if(FALSE){
## Tree of full-length SILVA sequences, indicating which mapped to GTDB
taxdata_seed <- taxdata_mothur %>%
  filter(!is.na(label)) %>%
    distinct( label, .keep_all=T) %>%
    select(c("label", "GTDBmapped", "Kingdom"))
seed.tree2 <- left_join(as_tibble(seed.tree1),taxdata_seed)

p <- 
  ggtree(as.treedata(seed.tree2), branch.length='none', layout='circular') + geom_tippoint(aes(color=Kingdom)) 

p

##Of seed db, what proportion are mapped?
seed.tree2$isTip <- isTip(seed.tree2, seed.tree2$node)
seed.tree2Tips <- seed.tree2 %>% filter(isTip==T)
}


if(FALSE){

t.v4 <- merge(as_tibble(seed.treev4) , taxdata, by.x="label", by.y="id") %>%
  filter(Kingdom=="Bacteria")





chry1 <- tree_subset(as.treedata(v4data), node=536, levels_back=0)
ggtree(chry1)+ geom_tippoint(aes(color=Genus, shape=isGenusNode)) + geom_nodepoint(aes(color=Genus, shape=isGenusNode))
chry1 <- chry1 %>% as_tibble

##Bacillaceae
p1 <- ggtree(BacNodes)+ geom_tippoint(aes(color=Family, shape=isFamilyNode, size=3)) +geom_nodepoint(aes(color=Family, shape=isFamilyNode, size=3)) 
p1 %>% collapse(node=709)

}


## Look within Phyla (defined by optimal threshold) 
## 0.5 isn't optimal for Phyla, but is close, and doesn't result in as many split phyla
tempData <- resultData$tax0.37
tempData$isPhylumNode <- is.na(tempData$Phylum[tempData$parent]) & !is.na(tempData$Phylum)
tempData <- tempData %>% mutate(Kingdom = na_if(Kingdom, ""),
                                Phylum = na_if(Phylum, ""),
                                Class = na_if(Class, ""),
                                Order = na_if(Order, ""),
                                Family = na_if(Family, ""),
                                Genus = na_if(Genus, ""),
                                Species = na_if(Species, ""))

## Pull out Phylum Nodes based on above cutoff
pNodes <- tempData %>% filter(isPhylumNode==T)
## Change NAs in InputData
inputData <- inputData %>% mutate(Kingdom = na_if(Kingdom, ""),
                                  Phylum = na_if(Phylum, ""),
                                  Class = na_if(Class, ""),
                                  Order = na_if(Order, ""),
                                  Family = na_if(Family, ""),
                                  Genus = na_if(Genus, ""),
                                  Species = na_if(Species, ""))

# Initialize outputs
subOutputScores <- data.frame()
subOutputNs <- data.frame()
subResultData <- list()

# For each Phylum, pull out relevent interior nodes beneath it
for(pNode in pNodes$node[which(!pNodes$isTip)]){
  phyCh <- offspring(tempData, pNode)
  if(sum(is.na(phyCh$Genus))==length(phyCh$Genus))
    next()
  
  #Save initial data for each cutoff, for later taxonomy assignments
  for(cut1 in cutoffs){
    subResultData[[paste0("tax",cut1)]] <- inputData
  }
  
  for (intNode in phyCh$node[which(phyCh$isTip==F)]){
    print(paste("Node:", intNode))
    ch <- offspring(inputData, intNode, tiponly=T)
    tre <- tree_subset(as.treedata(inputData), intNode, levels_back = 0)
    maxDist <- max(cophenetic(as.phylo(tre)))
    
    
    for(cut1 in cutoffs){
      if(maxDist < cut1) {
        if(length(table(ch[["Genus"]]))!=0)
          subResultData[[paste0("tax",cut1)]][intNode, "Genus"] <- names(table(ch[["Genus"]]))[which(table(ch[["Genus"]])==max(table(ch[["Genus"]])))][[1]]
      }
    }
  }
  
  for(cut1 in cutoffs){
    print(cut1)
    tempData2 <- subResultData[[paste0("tax",cut1)]]
    
    ### Define "Genus" as a clade defined by a Genus-named node with a non-Genus-named parent node, etc
    tempData2$isGenusNode <- is.na(tempData2$Genus[tempData2$parent]) & !is.na(tempData2$Genus)

    ## Calculate Genus Scores
    gNodes <- tempData2 %>% filter(isGenusNode==T & node %in% phyCh$node)
    ## First score to minimize: number of Genera represented multiple times
    ## Second score to minimize: number of non-Genus members clustered within a Genus
    genera <- offspring(tempData2, tempData2$node[which(tempData2$node %in% gNodes$node)], tiponly=T, self_include=T)
    if(is.null(genera[["Genus"]])){
      generaTables <- lapply(genera, FUN= function(x) table(x[["Genus"]]))
      generaNames <- tempData2 %>% filter(node %in% gNodes$node) %>% select(Genus)
      generaNames$score2count <-  unlist(lapply(generaTables, FUN=function(x) sum(x)-max(x)))
    } else {
      generaTables <- table(genera[["Genus"]])
      generaNames <- tempData2 %>% filter(node %in% gNodes$node) %>% select(Genus)
      generaNames$score2count <-  sum(generaTables)-max(generaTables)
      }
    
    
    generascore2s <- generaNames %>% group_by(Genus) %>% summarize(score2count = sum(score2count))
    generascore2s$score1count <- as.integer(table(gNodes$Genus)-1)
    gscore1 <- sum(generascore2s$score1count)
    gscore2 <- sum(generascore2s$score2count)
    gscoreSum <- gscore1 + gscore2
    subOutputScores[as.character(pNode), as.character(cut1)] <- gscoreSum/sum(unlist(lapply(genera,nrow)))
    subOutputNs[as.character(pNode), as.character(cut1)] <- sum(unlist(lapply(genera,nrow)))

  }
  
  nseqs <- nrow(inputData %>% filter(isTip==TRUE))
  
  
}
subOutputPlot <- 
  subOutputScores[which(rownames(subOutputScores) %in% rownames(subOutputNs)[which(subOutputNs[,1]>=10)]),]

nodes <- as.numeric(rownames(subOutputPlot))

subOutputPlot <- subOutputPlot[rev(order(pull(tempData, "Phylum")[nodes])),]
subOutputPlotNames <- pull(tempData, "Phylum")[nodes][rev(order(pull(tempData, "Phylum")[nodes]))]

subOutputPlotLong <- subOutputPlot %>% 
  tibble::rownames_to_column(var = "node") %>%
  pivot_longer( -node, names_to = "Threshold", values_to = "Error_Score") %>%
  group_by(node) %>%
  filter(Error_Score==min(Error_Score))

subOutputPlotLong$node <- as.integer(subOutputPlotLong$node)

subOutputPlotLong <- left_join(subOutputPlotLong, tempData, by="node") %>%
  arrange(desc(Phylum))
subOutputPlotLong$PhylumNode <- paste0(subOutputPlotLong$Phylum, subOutputPlotLong$node)

subOutputNs2 <- data.frame(node=as.integer(rownames(subOutputNs)), Nseqs=subOutputNs[,1])
subOutputPlotLong <- left_join(subOutputPlotLong, subOutputNs2, by="node")

p<-ggplot(subOutputPlotLong, aes(x=PhylumNode, y=as.numeric(Threshold), size=Nseqs)) + 
  geom_point() + ylim(0, 0.37) + 
  theme_gray(base_size = 14) +
  ylab("Best Genus Threshold") + xlab("Phylum Node") + geom_hline(yintercept=0.04)
p + coord_flip()

p<-ggplot(subOutputPlotLong, aes(x=Phylum, y=as.numeric(Threshold), size=Nseqs)) + 
  geom_point() + ylim(0, 0.37) + 
  ylab("Best Genus Threshold") + xlab("Phylum Node") + geom_hline(yintercept=0.04)
p + coord_flip()



test <- tree_subset(as.treedata(inputData), 6390, levels_back = 0)
ggtree(test) + geom_tippoint(aes(color=Genus))

test <- tree_subset(as.treedata(v4treedata), 9015, levels_back = 2) %>% as_tibble
ggtree(test) + geom_tippoint(aes(color=Genus))

save(resultData, file="~/proj/16s-region-checker/output/MappingBasedOnBinomialModel.RData")
#resultDataBin <- resultData
#load("~/proj/16s-region-checker/output/MappingBasedOnErrorScores.RData")
