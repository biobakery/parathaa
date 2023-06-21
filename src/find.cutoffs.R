#!/usr/bin/env Rscript

### This script finds optimal cutoffs for each taxonomic level for a given tree, and plots their error scores

require(docopt)
'Usage:
   find.cutoffs.R [-d <naming file> -o <output> -n <tree>]

Options:
   -d naming file [default: input/taxmap_slv_ssu_ref_138.1.txt]
   -o model output directory [default: output/testrun20230117]
   -n tree [default: /Users/mis696/proj/parathaa/output/20230109_SyntheticV4V5_nameHarmonizing/region_specific.tree]
   --wt1 over-split penalty weight [default: 1]
   --wt2 over-merge penalty weight [default: 1]

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
##library(ggimage)
library(ggplot2)
library(TDbook)
source("src/SILVA.species.editor.R")
source("src/calc.error.scores.R")

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

taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")

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
in.tree.data <- SILVA.species.editor(in.tree.data, Task="find_cutoffs")

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
  
}


# can we vectorize this loop?

## Go through the cutoffs
for(cut1 in cutoffs){
  print(cut1)
  tempData <- resultData[[paste0("tax",cut1)]]
  
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
