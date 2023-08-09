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
## load in the internal node stats we calculated in the previous job.
load(file.path(opts$o, "internal_node_stats.RData"))

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

pb = txtProgressBar(min = 0, max = length(internal_node_stats), initial = 0) 

for (i in 1:length(internal_node_stats)){
  setTxtProgressBar(pb,i)
  #progress bar?
  intNode <- which(inputData$isTip==F)[i]
  
  maxDist <- internal_node_stats[[i]][[1]]
  
  ch <- internal_node_stats[[i]][[2]]
  
  #print(paste("Node:", intNode))

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
}
close(pb)

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
