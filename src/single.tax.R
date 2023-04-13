## Assigns taxonomy for a single interior node at a single level, based on:
## intNode: interior node of a phylogenetic tree
## level: one of "Phylum", "Class", "Order", "Family", "Genus", "Species" at which to assign a taxonomic classification
## maxDist: maxiumum cophenetic distance among tips underneath intNode
## cutoff: optimal distance cutoffs for the taxonomic level
## nodeGroups: table of tip nodes under intNode at the given taxonomic level

single.tax <- function(intNode, level, maxDist, cutoff, nodeGroups, falseNegRate, acceptableProb, resultData){
  ## Initialize:
  results <- list()
  ## if the maximum distance among the child (tip) nodes is less than the cutoff for the taxonomic level,
  ## (and if there is at least one child node; aka the node is not a tip)
  ## then 
  if(maxDist < cutoff & length(nodeGroups)!=0){
    if(pbinom(sum(nodeGroups)-max(nodeGroups)-1,
              sum(nodeGroups),falseNegRate, lower.tail = F) > acceptableProb){
      
      # if that is correct we then assign the int node to that phylum
      resultData[["tax_bestcuts"]][intNode, level] <- names(nodeGroups)[which(nodeGroups==max(nodeGroups))][[1]]
    }
    # if it doesn't pass we assign multiple names to that node 
    # were names are base on the assignment that passes that error model.
    else {resultData[["tax_bestcuts"]][intNode, level] <- paste(
      names(which(nodeGroups >
                    qbinom(acceptableProb , sum(nodeGroups),falseNegRate, lower.tail = F))), collapse=";")}
  }
    return(results)
}
