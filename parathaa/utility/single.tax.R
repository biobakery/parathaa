## Assigns taxonomy for a single interior node at a single level, based on:
## intNode: interior node of a phylogenetic tree
## level: one of "Phylum", "Class", "Order", "Family", "Genus", "Species" at which to assign a taxonomic classification
## maxDist: maxiumum cophenetic distance among tips underneath intNode
## cutoff: optimal distance cutoffs for the taxonomic level
## nodeGroups: table of tip nodes under intNode at the given taxonomic level

single.tax <- function(intNode, level, maxDist, cutoff, nodeGroups, falseNegRate, acceptableProb){
  ## Initialize:
  result <- NA

  #if its below the cutoff and there is assigned taxonomy to the tips do the following:
  if(maxDist < cutoff & length(nodeGroups)!=0){
    #check if there is only 1 dominant taxon and if so just return that label
    if(pbinom(sum(nodeGroups)-max(nodeGroups)-1,
              sum(nodeGroups),falseNegRate, lower.tail = F) > acceptableProb){
      # if that is correct we then assign the int node to that phylum
      result <- names(nodeGroups)[which(nodeGroups==max(nodeGroups))][[1]]
    }
    # if it doesn't pass we assign multiple names to that node 
    # were names are base on the assignment that passes that error model.
    else {result <- paste(
      names(which(nodeGroups >
                    qbinom(acceptableProb , sum(nodeGroups),falseNegRate, lower.tail = F))), collapse=";")}
    
    if(length(result)==0){
      result <- names(nodeGroups)
    }
  #if its below the maxDist and the below nodes don't have taxonomy then we can return "Unclassified"  
  }else if(maxDist < cutoff & length(nodeGroups)==0){
    result <- "Unclassified"
  }
  return(result)
}
