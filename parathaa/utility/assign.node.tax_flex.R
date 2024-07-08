#!/usr/bin/env Rscript


require(docopt)
'Usage:
   assign.node.tax.R [-d <naming file> -o <output> -n <tree> --bError <error_rate> --bThreshold <threshold> --util1 <spec.editor_PATH> --util2 <single.tax.R PATH>]

Options:
   -d naming file [default: input/taxmap_slv_ssu_ref_138.1.txt]
   -o model output directory [default: output/20230406_testrun]
   -n tree [default: output/20230406_testrun/region_specific.tree]
   --bError binomial error rate [default: 0.05]
   --bThreshold binomial error rate [default: 0.20]
   --util1 PATH to SILVA.species.editor.dev.R [default: utility/SILVA.species.editor.dev.R]
   --util2 PATH to single.tax.R [default: utility/single.tax.r] 
   
 ]' -> doc

#To add as arguments: binom error params

opts <- docopt(doc)

# If we wanted we could add a loop that tries to install these packages within R?
# Although it might be nice to have this as a seperate R script 
# that runs on anadama workflow start up

library(ggtree)
library(treeio)
library(tidyr)
library(dplyr)
library(ape)
library(TDbook)
source(opts$util1)
source(opts$util2)


hierachry <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
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


## we need to update this to read in a flexible taxonomy file...

taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="", check.names = FALSE)

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


### so at this point we remove all the species labels that say unknown and leave them blank which is what is causing the bug in the first place..
### if we allow unkowns then it wouldn't be an issue of having unknowns i don't think... but what other issues might this generate
### on the back end?

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

#pb = txtProgressBar(min = 0, max = length(internal_node_stats), initial = 0, style=3) 

for (i in 1:length(internal_node_stats)){
  #setTxtProgressBar(pb,i)
  #progress bar?
  intNode <- which(inputData$isTip==F)[i]
  
  maxDist <- internal_node_stats[[i]][[1]]
  
  ch <- internal_node_stats[[i]][[2]]
  
  #print(paste("Node:", intNode))

  # cutoffs represent distances on the tree from tip to tip
  
  resultData[["tax_bestcuts"]][intNode, "maxDists"] <- maxDist
  

  
  lastlabel <- ""
  for(level in hierachry){
    cutoff <- cutoffs[level]
    
    #should we treat this in the same way so we commit to this genus so we only consider things under that genus
    #if so we do something like this:
    
    
    ## here check if there is an ambiguous label
    if(grepl(";", lastlabel)){
      choosen_taxon <- str_split(lastlabel, ";")
      choosen_labels <- c()
      for(taxon in choosen_taxon[[1]]){
        #remove unclassified from the previous label if it exists for search purposes within the tree
        taxon_search <- gsub(" Unclassified", "", taxon)
        #grab the previous level in taxonomy
        previous_level <- hierachry[which(hierachry==level)-1]
        #filter to only tips of the chosen taxon
        taxon_filt_tab <- ch %>% filter(!!sym(previous_level)==taxon_search)
        #apply binomial error model
        nodeGroups <- table(taxon_filt_tab[[level]])
        ## apply binomial error model code
        label <- single.tax(intNode=intNode,
                            level=level, 
                            maxDist=maxDist, 
                            cutoff=cutoff, 
                            nodeGroups=nodeGroups, 
                            falseNegRate=falseNegRate, 
                            acceptableProb=acceptableProb)
        ##need to think about how to deal with unclassified...
        if(!is.na(label)){
         if(label=="Unclassified"){
           label <- paste(taxon_search, label, sep=" ")
         } 
        }
        choosen_labels <- c(choosen_labels, label)
      }
      ## if it fails for all then the final label should be NA
      if(length(which(is.na(choosen_labels)))==length(choosen_labels)){
        final_lab <- NA
      }else{
        final_lab <- paste0(choosen_labels, collapse=";")
      }
      
      resultData[["tax_bestcuts"]][intNode, level] <- final_lab
      lastlabel <- final_lab
    }else{
      #then we need to use the information from the last accepted assignment
      nodeGroups <- table(ch[[level]]) ## this is done here rather than inside the function single.tax 
      ## because we only want to calculate ch once for a node (it is time-consuming)
      
      label  <- single.tax(intNode=intNode,
                                                                 level=level, 
                                                                 maxDist=maxDist, 
                                                                 cutoff=cutoff, 
                                                                 nodeGroups=nodeGroups, 
                                                                 falseNegRate=falseNegRate, 
                                                                 acceptableProb=acceptableProb)
      if(!is.na(label)){
        if(label=="Unclassified"){
          label <- paste0(gsub(" Unclassified", "", lastlabel)," Unclassified")
        }
      }
      resultData[["tax_bestcuts"]][intNode, level] <- label
      #keep track of the last label we assigned
      lastlabel <- resultData[["tax_bestcuts"]][intNode, level]
    }
    
  }
}

#check if its assigned at X level
# check if its parent is not
# if so call it a XNode
resultData[["tax_bestcuts"]]$isKingdomNode <- is.na(resultData[["tax_bestcuts"]]$Kingdom[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Kingdom)
resultData[["tax_bestcuts"]]$isPhylumNode <- is.na(resultData[["tax_bestcuts"]]$Phylum[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Phylum)
resultData[["tax_bestcuts"]]$isClassNode <- is.na(resultData[["tax_bestcuts"]]$Class[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Class)
resultData[["tax_bestcuts"]]$isOrderNode <- is.na(resultData[["tax_bestcuts"]]$Order[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Order)
resultData[["tax_bestcuts"]]$isFamilyNode <- is.na(resultData[["tax_bestcuts"]]$Family[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Family)
resultData[["tax_bestcuts"]]$isGenusNode <- is.na(resultData[["tax_bestcuts"]]$Genus[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Genus)
resultData[["tax_bestcuts"]]$isSpeciesNode <- is.na(resultData[["tax_bestcuts"]]$Species[resultData[["tax_bestcuts"]]$parent]) & !is.na(resultData[["tax_bestcuts"]]$Species)

## also need to propagate unclassified labels to tips but we do this after internal label assignment
## this is because we don't want to consider unclassified as a label at the assigned level...

## fix taxonomy within ch but this adds a ton of run time... (but we don't really want to fix them for distance testing purposes...)
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Kingdom = ifelse(is.na(Kingdom) & isTip, "Unclassified", Kingdom),
                    Phylum = ifelse(is.na(Phylum) & isTip, paste0(Kingdom, "Unclassified", sep=" "), Phylum))

resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Class = ifelse(is.na(Class) & isTip, paste(Phylum,"Unclassified", sep=" "), Class))

#Fix Order
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Order = ifelse(grepl("Unclassified", Class) & isTip, Class, Order))
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Order = ifelse(is.na(Order) & isTip, paste(Class,"Unclassified", sep=" "), Order))

#Fix Family
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Family = ifelse(grepl("Unclassified", Order) & isTip, Order, Family))
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Family = ifelse(is.na(Family) & isTip, paste(Order,"Unclassified", sep=" "), Family))

#Fix Genus
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Genus = ifelse(grepl("Unclassified", Family) & isTip, Family, Genus))
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Genus = ifelse(is.na(Genus) & isTip, paste(Family,"Unclassified", sep=" "), Genus))

#Fix Species
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Species = ifelse(grepl("Unclassified", Genus) & isTip, Genus, Species))
resultData$tax_bestcuts <- resultData$tax_bestcuts %>% mutate(Species = ifelse(is.na(Species) & isTip, paste(Genus,"Unclassified", sep=" "), Species))

## Save output
save(resultData, file = file.path(opts$o, "resultTree_bestThresholds.RData"))




