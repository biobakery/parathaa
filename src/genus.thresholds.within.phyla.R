## The script below was taken from make.taxonomy.trees.R and was an early script to plot
## "best" genus-defining thresholds within each phylum, to assess variability across the ToL
## This script needs work, but I don't want to fully delete it because I plan to adapt it for
## a figure in the manuscript. Doesn't belong in the Parathaa pipeline, though.


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
      generaNames <- tempData2 %>% filter(node %in% gNodes$node) %>% dplyr::select(Genus)
      generaNames$score2count <-  unlist(lapply(generaTables, FUN=function(x) sum(x)-max(x)))
    } else {
      generaTables <- table(genera[["Genus"]])
      generaNames <- tempData2 %>% filter(node %in% gNodes$node) %>% dplyr::select(Genus)
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


save(resultData, file="~/proj/16s-region-checker/output/MappingBasedOnBinomialModel.RData")
#resultDataBin <- resultData
#load("~/proj/16s-region-checker/output/MappingBasedOnErrorScores.RData")

