#!/usr/bin/env Rscript

# General comment:
# It looks like this script has two different uses? (maybe part of legacy analysis?)
# If so I would suggest keep that in a seperate script and cleaning this so that it only
# has the parts that run in the anadama workflow...
require(docopt)
'Usage:
   make.taxonomy.trees.R [-d <naming file> -o <output> -t <task> -n <tree>]

Options:
   -d naming file [default: input/taxmap_slv_ssu_ref_138.1.txt]
   -o model output directory [default: output/testrun20230117]
   -t task: assign_Tax or find_cutoffs [default: find_cutoffs]
   -n tree [default: /Users/mis696/proj/parathaa/output/20230109_SyntheticV4V5_nameHarmonizing/region_specific.tree]

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
library(ggimage)
library(TDbook)
library(ggplot2)
source("src/SILVA.species.editor.R")

## Set task
Task <- opts$t #"assign_Tax" "find_cutoffs"

## Bring in taxonomy file
inFileTaxdata <- opts$d

## Tree made from database trimmed to region
in.tree <- read.newick(opts$n)
in.tree.data <- as_tibble(in.tree)
suppressWarnings({
in.tree.data <- in.tree.data %>%
  separate(col=label, into=c("primaryAccession", "arbID"), remove=F)
})
## The above gives a warning because only the tip nodes have primary accessions and arbIDs.
## I suppressed it so that the warning message doesn't reach the user.


## Taxonomy of database from SILVA (primary accession, start, stop)

# I wonder if using a data.table/tibbly/big.data (pckge) would increase performance here?
# (below three lines run fairly slowly...)

taxdata <- data.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")

# the IDs that come with the file specific ID to the sequence but they also have start and stop for the sequence!
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)

suppressWarnings({
taxdata <- taxdata %>%
  separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  dplyr::rename(Species = organism_name) %>%
  filter(Kingdom=="Bacteria") 
})

#in.tree.data is from the FAstTree... and taxdata is the silva taxonomy database
in.tree.data <- left_join(in.tree.data, taxdata, by="primaryAccession", multiple="all") %>%
  distinct(node, .keep_all=T) 

#some functions won't work unless its a specific type of class of data
class(in.tree.data) <- c("tbl_tree", class(in.tree.data))

#Check if each node is a tip or not!
in.tree.data$isTip <- isTip(in.tree.data, in.tree.data$node)

## Fix species names

# there are a lot of species names that are inconsistent and needed to be cleaned up!
in.tree.data <- SILVA.species.editor(in.tree.data, Task=Task)

in.tree.data <- in.tree.data %>% mutate(Kingdom = na_if(Kingdom, ""),
                                Phylum = na_if(Phylum, ""),
                                Class = na_if(Class, ""),
                                Order = na_if(Order, ""),
                                Family = na_if(Family, ""),
                                Genus = na_if(Genus, ""),
                                Species = na_if(Species, ""))  


# parameters from the binomModel
falseNegRate <- 0.05
acceptableProb <- 0.20


#Broad range of cutoffs:
if(Task=="find_cutoffs"){
  #looks at a series of cutoffs for the best "threshold"
  #might be better to look minima and choose that and search around it rather than looking at all these...
  cutoffs<-c(seq(0.001, 0.009, by=0.001), seq(0.01, 0.5, by=0.01), seq(0.55, 0.9, 0.05))
}
#"Best" cutoffs loaded from previous step:
if(Task=="assign_Tax"){
  load(file.path(opts$o, "optimal_scores.RData"))
  bestThresh <- plotData2 %>% group_by(Level) %>% summarise(minThreshold = mean(minThreshold))
  cutoffs <- bestThresh$minThreshold
  names(cutoffs) <- bestThresh$Level
  #cutoffs<-c("Species"=0.003, "Genus"=0.06, "Family"=0.13, "Order"=0.21,  "Class"=0.36, "Phylum"=0.46)
}

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
if(Task=="find_cutoffs"){
  for(cut1 in cutoffs){
    ## Note from Meg: This step is not currently memory-efficient, it stores a lot of the same data many times. 
    resultData[[paste0("tax",cut1)]] <- inputData
  }
}
if(Task=="assign_Tax")
  resultData[["tax_bestcuts"]] <- inputData
startTime <- Sys.time()

### This is a very complicated chunk of code and I'm guessing it is the "meat"
### of this script. I wonder if we can break up the tasks a bit more into functions
### That would make the code a bit more readable 
### We can then also do unit tests on each function to make sure they all work as expected!


## In general it would be helpful to comment on each few lines below to describe what they are doing
## I did my best to write my best guess as to what each line of code does without deep diving


# profile the tool to figure out where its using the most memory and see whats going on there


# Also a general comment! For loops and very slow in R I wonder if we can make use of other
# If we cannot vectorize this loop (using apply) maybe we can think about allow multi core
# loop using the "foreach" package?

## go through each node that is not a tip
for (intNode in which(inputData$isTip==F)){
  #progress bar?
  print(paste("Node:", intNode))
  #grab tips
  ch <- offspring(inputData, intNode, tiponly=T)
  #get the sub tree
  tre <- tree_subset(as.treedata(inputData), intNode, levels_back = 0)
  # grab the max cophenetic distance of the tips within the sub tree?
  # I'm guessing this is what is used to determine the cut-off of if we can be sure its a specific taxa?
  
  maxDist <- max(cophenetic(as.phylo(tre)))
  
  
  # cutoffs represent distances on the tree from tip to tip
  
  # do a correlation plot from distances to ANI would be nice!
  
  if(Task=="find_cutoffs"){
    for(cut1 in cutoffs){
      # go through at each level and check if the distance is less than the cut-off
      # help determine which cut-off performs the best in the dataset?
      # unclear how the "best" cut-off is measured at this point.
      # if maxDist is greater than the cutoff what happens?
          # for each cut point you want to calculate the error
          # two parts over clumping and over splitting
      # 
      if(maxDist < cut1) {
        #as long as the children have been assigned a kingdom (what below if statement does)
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
    # in resultData many of maxdist as NA why is that?
    # just grab this warning if we can
    
  } else if(Task=="assign_Tax"){
    resultData[["tax_bestcuts"]][intNode, "maxDists"] <- maxDist
    
    #check that atleast one child node has been assigned at this taxonomy level
    if(length(table(ch[["Kingdom"]]))!=0)
      
      #maximum name which doesn't really matter with Kingdom
      resultData[["tax_bestcuts"]][intNode, "Kingdom"] <- names(table(ch[["Kingdom"]]))[which(table(ch[["Kingdom"]])==max(table(ch[["Kingdom"]])))][[1]]
    #it seems the below line of code is run in the same manner but on each level of taxonomy
    #this would be a good place to convert it into a function (with good documentation)
    #and use the level of taxonomy as part of the input for that function!
    
    # the distance has to be lower than the optimal cutoff and it has children node that were assigned phylum
    # if its above cutoff we leave it as "unassigned/unclassified/however you want to call it"
    if(maxDist < cutoffs["Phylum"] & length(table(ch[["Phylum"]]))!=0){
      #i'm guessing this is checking whether this taxonomy level meets some
      # specification within its binomial distribution
      # q is calculated in the first line,
      # size is than the sum of the total number of children within that phylum?
      # falsenegrate (constant at 0.05) is used as the prob parameter?
      # then we check if that is greater than the accetableProb 
      # not sure how that is defined and why its a contanst of 0.2
      
      
      # okay so quantiles = the number of that phylum - the majority number
      # if the number of nodes that are not the majority phylum are larger than 
      # the error threhold would we accept that error rate
      
      # check if all of the children are one phylum? (might speed it up??)
      if(pbinom(sum(table(ch[["Phylum"]]))-max(table(ch[["Phylum"]]))-1,
                sum(table(ch[["Phylum"]])),falseNegRate, lower.tail = F) > acceptableProb){
        
        # if that is correct we then assign the int node to that phylum
        resultData[["tax_bestcuts"]][intNode, "Phylum"] <- names(table(ch[["Phylum"]]))[which(table(ch[["Phylum"]])==max(table(ch[["Phylum"]])))][[1]]
      } 
        # if it doesn't pass we assign multiple names to that node 
        # were names are base on the assignment that passes that error model.
      else {resultData[["tax_bestcuts"]][intNode, "Phylum"] <- paste(
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

  } 
    
}

if(Task=="assign_Tax"){
  ### Define "Genus" as a clade defined by a Genus-named node with a non-Genus-named parent node, etc
  
  #check if its assigned at X level
  # check if its parent is not
  # if so call it a XNode
  resultData[["tax_bestcuts"]]$isPhylumNode <- is.na(resultData[["tax_bestcuts"]]$Phylum[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Phylum)
  resultData[["tax_bestcuts"]]$isClassNode <- is.na(resultData[["tax_bestcuts"]]$Class[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Class)
  resultData[["tax_bestcuts"]]$isOrderNode <- is.na(resultData[["tax_bestcuts"]]$Order[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Order)
  resultData[["tax_bestcuts"]]$isFamilyNode <- is.na(resultData[["tax_bestcuts"]]$Family[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Family)
  resultData[["tax_bestcuts"]]$isGenusNode <- is.na(resultData[["tax_bestcuts"]]$Genus[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Genus)
  resultData[["tax_bestcuts"]]$isSpeciesNode <- is.na(resultData[["tax_bestcuts"]]$Species[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Species)
}

elapsedTime <- Sys.time() - startTime
print(paste("Time elapsed:", elapsedTime, "minutes"))

# identitation for the below code is off making is hard to tell where this loop
# starts and ends


if(Task=="find_cutoffs"){
  # can we vectorize this loop?
  
  #go through the cutoffs
  for(cut1 in cutoffs){
    print(cut1)
    tempData <- resultData[[paste0("tax",cut1)]]
    # taking the list item
    # for each int node its assigned the majority class at each level
    
    ### Define "Genus" as a clade defined by a Genus-named node with a non-Genus-named parent node, etc
    
    # we are trying to find intNodes that we can consider everything under it to be of the same class
    # we are looking for intnodes where its parent has children from multiple taxa 
    # we then define those as "genus node etc."
    tempData$isKingdomNode <- is.na(tempData$Kingdom[tempData$parent]) & !is.na(tempData$Kingdom)
    tempData$isPhylumNode <- is.na(tempData$Phylum[tempData$parent]) & !is.na(tempData$Phylum)
    tempData$isClassNode <- is.na(tempData$Class[tempData$parent]) & !is.na(tempData$Class)
    tempData$isOrderNode <- is.na(tempData$Order[tempData$parent]) & !is.na(tempData$Order)
    tempData$isFamilyNode <- is.na(tempData$Family[tempData$parent]) & !is.na(tempData$Family)
    tempData$isGenusNode <- is.na(tempData$Genus[tempData$parent]) & !is.na(tempData$Genus)
    tempData$isSpeciesNode <- is.na(tempData$Species[tempData$parent]) & !is.na(tempData$Species)
    
    
    ## below looks like another nice code block that could probably be turned into a function
    ## since its called at multiple taxonomic levels...
    ## Calculate Phylum Scores
    pNodes <- tempData %>% filter(isPhylumNode==T)
    ## First score to minimize: number of Phylum represented multiple times
    ## Second score to minimize: number of non-Phylum members clustered within a Phylum
    
    # grab all the Phylum nodes from Tempdata
    phyla <- offspring(tempData, tempData$node[which(tempData$isPhylumNode==T)], tiponly=T, self_include=T)
    
    # for each phylum node we go through and calculate the number of phylum in that phyla node
    phylumTables <- lapply(phyla, FUN= function(x) table(x[["Phylum"]]))
    
    # we go and select the nodes that are counted above
    phylumNames <- tempData %>% filter(node %in% names(phylumTables)) %>% dplyr::select(Phylum, node)
    
    #taking total number of phyla within the phylum node and - max which is the majority class and name of that node
    phylumNames$score2count <-  unlist(lapply(phylumTables, FUN=function(x) sum(x)-max(x)))
    
    #summarize this score across all phyla
    phylumscore2s <- phylumNames %>% group_by(Phylum) %>% summarize(score2count = sum(score2count))
    
    # number of phyla reprented multiple times
    # goes through the name of the phylum and get see if there are multiple and scores it (first score above)
    phylumscore2s$score1count <- as.integer(table(pNodes$Phylum)-1)
    pscore1 <- sum(phylumscore2s$score1count)
    pscore2 <- sum(phylumscore2s$score2count)
    pscoreSum <- pscore1 + pscore2
    phylumscore2s$scoreSum <- phylumscore2s$score1count + phylumscore2s$score2count
    
    outputCounts[["Phylum"]][[as.character(cut1)]] <- phylumscore2s
    outputScores[["Phylum"]][as.character(cut1)] <- pscoreSum
    #outputNs[["Phylum"]][as.character(cut1)] <- (lapply(phyla,nrow))
    
    ## Calculate Class Scores
    cNodes <- tempData %>% filter(isClassNode==T)
    ## First score to minimize: number of Class represented multiple times
    ## Second score to minimize: number of non-Class members clustered within a Class
    classes <- offspring(tempData, tempData$node[which(tempData$isClassNode==T)], tiponly=T, self_include=T)
    classTables <- lapply(classes, FUN= function(x) table(x[["Class"]]))
    classNames <- tempData %>% filter(node %in% names(classTables)) %>% dplyr::select(Class, node)
    classNames$score2count <-  unlist(lapply(classTables, FUN=function(x) sum(x)-max(x)))
    classscore2s <- classNames %>% group_by(Class) %>% summarize(score2count = sum(score2count))
    classscore2s$score1count <- as.integer(table(cNodes$Class)-1)
    cscore1 <- sum(classscore2s$score1count)
    cscore2 <- sum(classscore2s$score2count)
    cscoreSum <- cscore1 + cscore2
    classscore2s$scoreSum <- classscore2s$score1count + classscore2s$score2count
    
    outputCounts[["Class"]][[as.character(cut1)]] <- classscore2s
    outputScores[["Class"]][as.character(cut1)] <- cscoreSum
    #outputNs[["Class"]][as.character(cut1)] <- (lapply(classes,nrow))
    
    ## Calculate Order Scores
    oNodes <- tempData %>% filter(isOrderNode==T)
    ## First score to minimize: number of Orders represented multiple times
    ## Second score to minimize: number of non-Order members clustered within a Order
    orders <- offspring(tempData, tempData$node[which(tempData$isOrderNode==T)], tiponly=T, self_include=T)
    orderTables <- lapply(orders, FUN= function(x) table(x[["Order"]]))
    orderNames <- tempData %>% filter(node %in% names(orderTables)) %>% dplyr::select(Order, node)
    orderNames$score2count <-  unlist(lapply(orderTables, FUN=function(x) sum(x)-max(x)))
    orderscore2s <- orderNames %>% group_by(Order) %>% summarize(score2count = sum(score2count))
    orderscore2s$score1count <- as.integer(table(oNodes$Order)-1)
    oscore1 <- sum(orderscore2s$score1count) 
    oscore2 <- sum(orderscore2s$score2count)
    oscoreSum <- oscore1 + oscore2
    orderscore2s$scoreSum <- orderscore2s$score1count + orderscore2s$score2count
    
    outputCounts[["Order"]][[as.character(cut1)]] <- orderscore2s
    outputScores[["Order"]][as.character(cut1)] <- oscoreSum
    #outputNs[["Order"]][as.character(cut1)] <- (lapply(orders,nrow))
    
    ## Calculate Family Scores
    fNodes <- tempData %>% filter(isFamilyNode==T)
    ## First score to minimize: number of Family represented multiple times
    ## Second score to minimize: number of non-Family members clustered within a Family
    families <- offspring(tempData, tempData$node[which(tempData$isFamilyNode==T)], tiponly=T, self_include=T)
    familyTables <- lapply(families, FUN= function(x) table(x[["Family"]]))
    familyNames <- tempData %>% filter(node %in% names(familyTables)) %>% dplyr::select(Family, node)
    familyNames$score2count <- unlist(lapply(familyTables, FUN=function(x) sum(x)-max(x)))
    familyscore2s <- familyNames %>% group_by(Family) %>% summarize(score2count = sum(score2count))
    familyscore2s$score1count <- as.integer(table(fNodes$Family)-1)
    fscore1 <- sum(familyscore2s$score1count)
    fscore2 <- sum(familyscore2s$score2count)
    fscoreSum <- fscore1 + fscore2
    familyscore2s$scoreSum <- familyscore2s$score1count + familyscore2s$score2count
    
    outputCounts[["Family"]][[as.character(cut1)]] <- familyscore2s
    outputScores[["Family"]][as.character(cut1)] <- fscoreSum
    #outputNs[["Family"]][as.character(cut1)] <- (lapply(families,nrow))
    
    
    ## Calculate Genus Scores
    gNodes <- tempData %>% filter(isGenusNode==T)
    ## First score to minimize: number of Genera represented multiple times
    ## Second score to minimize: number of non-Genus members clustered within a Genus
    genera <- offspring(tempData, tempData$node[which(tempData$isGenusNode==T)], tiponly=T, self_include=T)
    generaTables <- lapply(genera, FUN= function(x) table(x[["Genus"]]))
    generaNames <- tempData %>% filter(node %in% names(generaTables)) %>% dplyr::select(Genus, node)
    generaNames$score2count <-  unlist(lapply(generaTables, FUN=function(x) sum(x)-max(x)))
    generascore2s <- generaNames %>% group_by(Genus) %>% summarize(score2count = sum(score2count))
    generascore2s$score1count <- as.integer(table(gNodes$Genus)-1)
    gscore1 <- sum(generascore2s$score1count)
    gscore2 <- sum(generascore2s$score2count)
    gscoreSum <- gscore1 + gscore2
    generascore2s$scoreSum <- generascore2s$score1count + generascore2s$score2count
    
    outputCounts[["Genus"]][[as.character(cut1)]] <- generascore2s
    outputScores[["Genus"]][as.character(cut1)] <- gscoreSum
    #outputNs[["Genus"]][as.character(cut1)] <- (lapply(genera,nrow))
    
    ## Calculate Species Scores
    sNodes <- tempData %>% filter(isSpeciesNode==T)
    ## First score to minimize: number of Species represented multiple times
    ## Second score to minimize: number of non-Species members clustered within a Species
    species <- offspring(tempData, tempData$node[which(tempData$isSpeciesNode==T)], tiponly=T, self_include=T)
    speciesTables <- lapply(species, FUN= function(x) table(x[["Species"]]))
    speciesNames <- tempData %>% filter(node %in% names(speciesTables)) %>% dplyr::select(Species, node)
    speciesNames$score2count <-  unlist(lapply(speciesTables, FUN=function(x) sum(x)-max(x)))
    speciesNames$nSeqs <- unlist(lapply(species, nrow))
    speciesscore2s <- speciesNames %>% group_by(Species) %>% summarize(score2count = sum(score2count))
    speciesscore2s$score1count <- as.integer(table(sNodes$Species)-1)
    sscore1 <- sum(speciesscore2s$score1count)
    sscore2 <- sum(speciesscore2s$score2count)
    sscoreSum <- sscore1 + sscore2
    speciesscore2s$scoreSum <- speciesscore2s$score1count + speciesscore2s$score2count
    outputCounts[["Species"]][[as.character(cut1)]] <- speciesscore2s
    outputScores[["Species"]][as.character(cut1)] <- sscoreSum
    #outputNs[["Species"]][as.character(cut1)] <- (lapply(species,nrow))
    
    
  }
}
  
if(Task=="find_cutoffs")
  save(resultData, file = file.path(opts$o, "resultTree_allThresholds.RData"))
if(Task=="assign_Tax")
  save(resultData, file = file.path(opts$o, "resultTree_bestThresholds.RData"))

print(Sys.time())


## Create Plot of Error Scores for each threshold at each level
if(Task=="find_cutoffs"){
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
  
  save(plotData2, file = file.path(opts$o, "optimal_scores.RData"))
  ggsave(filename = file.path(opts$o, "optimal_scores.png"), height = 4, width=5, units = "in")
}
