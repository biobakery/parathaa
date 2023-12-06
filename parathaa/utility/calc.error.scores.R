# This function calculates error scores for a given taxonomic level
## input is a tree data object and a taxonomic level, output is a list with 
## 1) data frame of error scores for each group at the taxonomic level, and 
## 2) sum of error scores across all groups at the taxonomic level

calc.error.scores <- function(treeData, level, wt1=1, wt2=1){
  #wt1 is the weight given to penalizing over-splitting
  #wt2 is the weight given to penalizing over-clumping
  
  results <- list()
  
  levelDefiningVar <- paste0("is", level, "Node")
  
  Nodes <- treeData %>% filter(!!as.symbol(levelDefiningVar)==T)
  
  #in some cases the cut-off might be too high and no nodes are assigned at this level
  #in this case we will return NA
  if(dim(Nodes)[1]==0){
    return(NULL)
  }
  ## First score to minimize: number of Phylum represented multiple times
  ## Second score to minimize: number of non-Phylum members clustered within a Phylum
  
  phylo_treeData <- as.phylo(treeData)
  groupsAtLevel <- list()
  isLeveltips <- treeData$node[which(treeData[[levelDefiningVar]]==T & treeData$isTip==T)]
  isLeveltips
  if(length(isLeveltips) != 0){
    for(i in 1:length(isLeveltips)){
      groupsAtLevel[[as.character(isLeveltips[i])]] <- treeData[which(treeData$node==isLeveltips[i]),]
    }
  }
  isLevelnode <- treeData$node[which(treeData[[levelDefiningVar]]==T & treeData$isTip==F)]
  if(length(isLevelnode) != 0){
    subNode_trees <- get_subtrees_at_nodes(phylo_treeData, isLevelnode-length(phylo_treeData$tip.label))
    for(i in 1:length(isLevelnode)){
      tips_to_grab <- subNode_trees$new2old_tips[[i]]
      groupsAtLevel[[as.character(isLevelnode[i])]] <- treeData[which(treeData$node %in% tips_to_grab),] 
    }
  }

  # for each group node we go through and calculate the number of names under that group node 
  # (aka for a clade that should be a genus based on the optimal distance threshold, how many different genera are under it)
  levelTables <- lapply(groupsAtLevel, FUN= function(x) table(x[[level]]))
  
  # we go and select the nodes that are counted above
  groupNames <- treeData %>% filter(node %in% names(levelTables)) %>% dplyr::select(!!as.symbol(level), node)
  
  #taking total number of names within the grouping node and - max which is the majority class and name of that node
  # ("over-clumping")
  groupNames$score2count <-  unlist(lapply(levelTables, FUN=function(x) sum(x)-max(x)))
  
  #sum this score across all groups with a given name (gives an error score for each name at a level)
  levelscore2s <- groupNames %>% group_by(!!as.symbol(level)) %>% summarize(score2count = sum(score2count))
  
  # number of names represented multiple times ("over-splitting")
  # goes through the name of the phylum and get see if there are multiple and scores it (first score above)
  levelscore2s$score1count <- as.integer(table(Nodes[[level]])-1)
  score1 <- sum(levelscore2s$score1count)
  score2 <- sum(levelscore2s$score2count)
  scoreSum <- ((wt1*score1) + (wt2*score2))/((wt1+wt2)/2)
  levelscore2s$scoreSum <- (wt1*levelscore2s$score1count + wt2*levelscore2s$score2count)/((wt1+wt2)/2)
  
  results[["counts"]] <- levelscore2s
  results[["scores"]] <- scoreSum
  
  return(results)
  
}
