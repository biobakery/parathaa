#!/usr/bin/env Rscript


require(docopt)
'Usage:
   assign.node.tax.R [-d <naming file> -o <output> -n <tree> --bError <error_rate> --bThreshold <threshold>]

Options:
   -d naming file [default: input/taxmap_slv_ssu_ref_138.1.txt]
   -o model output directory [default: output/20230406_testrun]
   -n tree [default: output/20230406_testrun/region_specific.tree]
   --bError binomial error rate [default: 0.05]
   --bThreshold binomial error rate [default: 1]
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
library(TDbook)
source("src/SILVA.species.editor.R")
source("src/single.tax.R")


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


## we need to update this to read in a flexible taxonomy file...

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
falseNegRate <- as.numeric(opts$bError)
acceptableProb <- as.numeric(opts$bThreshold)

## Get optimal scores from find.cutoffs.R output
load(file.path(opts$o, "optimal_scores.RData"))
bestThresh <- plotData2 %>% group_by(Level) %>% summarise(minThreshold = mean(minThreshold))
cutoffs <- bestThresh$minThreshold
names(cutoffs) <- bestThresh$Level
#cutoffs<-c("Species"=0.003, "Genus"=0.06, "Family"=0.13, "Order"=0.21,  "Class"=0.36, "Phylum"=0.46)


#Initialize
outputScores <- list()
outputNs <- list()
outputCounts <- list()
resultData <- list()
inputData <- in.tree.data
inputData$maxDists <- NA


print(acceptableProb)
resultData[["tax_bestcuts"]] <- inputData


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
  
  
  resultData[["tax_bestcuts"]][intNode, "maxDists"] <- maxDist
  
  #check that atleast one child node has been assigned at this taxonomy level
  if(length(table(ch[["Kingdom"]]))!=0){
    #maximum name which doesn't really matter with Kingdom
    resultData[["tax_bestcuts"]][intNode, "Kingdom"] <- names(table(ch[["Kingdom"]]))[which(table(ch[["Kingdom"]])==max(table(ch[["Kingdom"]])))][[1]]
  }
  #it seems the below line of code is run in the same manner but on each level of taxonomy
  #this would be a good place to convert it into a function (with good documentation)
  #and use the level of taxonomy as part of the input for that function!
  
  for(level in c("Phylum", "Class", "Order", "Family", "Genus", "Species")){
    cutoff <- cutoffs[level]
    nodeGroups <- table(ch[[level]]) ## this is done here rather than inside the function single.tax 
    ## because we only want to calculate ch once for a node (it is time-consuming)
  
    resultData[["tax_bestcuts"]][intNode, level] <- single.tax(intNode=intNode,
                                                               level=level, 
                                                               maxDist=maxDist, 
                                                               cutoff=cutoff, 
                                                               nodeGroups=nodeGroups, 
                                                               falseNegRate=falseNegRate, 
                                                               acceptableProb=acceptableProb
                                                               )
  
  }
  
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


## Save output
save(resultData, file = file.path(opts$o, "resultTree_bestThresholds.RData"))
