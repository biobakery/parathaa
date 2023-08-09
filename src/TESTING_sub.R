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

## we can write this loop in two parts...
## we can first write a parallel loop that will return the maxDist and the subtree tips in a list with two elements
## we can then take each one of those


start_time1 <- Sys.time()
results <- foreach(i=1:length(cutoffs), .combine = c) %dopar% {
  cut1 <- cutoffs[[i]]
  for(intNode in which(inputData$isTip==F)){
    intNode <- which(inputData$isTip==F)[i]
    print(paste("Node:", intNode))
    
    tmp_node <- intNode - length(which(inputData$isTip==T))
    sub_tree <- castor::get_subtree_at_node(treeio::as.phylo(inputData), tmp_node)
    maxDist <- castor::find_farthest_tip_pair(sub_tree$subtree)$distance
    
    tmp_tips <- sub_tree$new2old_tip
    ch <- inputData[which(inputData$node %in% tmp_tips),]
  
    tmp_rez <- list()
    tmp_rez[[paste0("tax",cut1)]] <- inputData
  
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
      if(length(table(ch[["Species"]]))!=0 & cut1 < 0.05)
        tmp_rez[[paste0("tax",cut1)]][intNode, "Species"] <- names(table(ch[["Species"]]))[which(table(ch[["Species"]])==max(table(ch[["Species"]])))][[1]]
    }
  }
  return(tmp_rez)
}
end_time1 <- Sys.time()


start_time2 <- Sys.time()
internal_node_stats <- foreach(i=1:length(which(inputData$isTip==F))) %do% {
  tmp_list <- list()
  
  intNode <- which(inputData$isTip==F)[i]
  print(paste("Node:", intNode))
  
  tmp_node <- intNode - length(which(inputData$isTip==T))
  sub_tree <- castor::get_subtree_at_node(treeio::as.phylo(inputData), tmp_node)
  maxDist <- castor::find_farthest_tip_pair(sub_tree$subtree)$distance
  
  tmp_tips <- sub_tree$new2old_tip
  ch <- inputData[which(inputData$node %in% tmp_tips),]
  
  tmp_rez <- list()
  tmp_rez[[paste0("tax",cut1)]] <- inputData
  
  tmp_list[["maxDist"]] <- maxDist
  tmp_list[["tips"]] <- ch
  
  return(tmp_list)
}
end_time2 <- Sys.time()
#took 21 seconds

start_time3 <- Sys.time()
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
end_time3 <- Sys.time()
