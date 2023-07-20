get_Max_cophenetic <- function(tree, node){
  max_dist <- 0
  
  tmp_tib <- as_tibble(tree)
  tmp_tib$isTip <- isTip(tmp_tib, tmp_tib$node)
  
  for(i in test_tib$node[which(test_tib$is_tip)]){
    message(i)
    tmp1 <- nodepath(as.phylo(tre), node, i)
    list_of_edges <- c()
    for(i in 2:length(tmp1)){
      list_of_edges <- c(list_of_edges, which(test$edge[,2] == tmp1[i]))
    }
    tmp_dist <- sum(test$edge.length[list_of_edges])
    if(tmp_dist > max_dist){
      max_dist <- tmp_dist
      message(max_dist)
    }

  }
  return(max_dist)
}
#there might be some kind of heuristic you can write where we go through and only sum them if there are many nodes?
#hmmm not really sure what to do about that
#at this point it would take 72 hours to calculate it for one internal 65 node...
# we then have 30,365 internal nodes to test...
#now each internal node will have a lower number of paths but not a lot less...
#

length(which(!test_tib$is_tip))


test_max_res <- get_Max_cophenetic(tre, intNode)

## need to figure out where to go from here but should go eat first...


#we need to get a list of tips
test_tib$is_tip <- isTip(test_tib, test_tib$node)
test_tib$node[which(test_tib$is_tip)]


test <- as.phylo(tre)
tmp1 <- nodepath(as.phylo(tre), 137647, 100)

#contains the path of nodes to the root
tmp1

#list of edges we need
list_of_edges <- c()
for(i in 2:length(tmp1)){
  list_of_edges <- c(list_of_edges, which(test$edge[,2] == tmp1[i]))
}
tmp1
sum(test$edge.length[list_of_edges])
#distance is 0.02106

which(test$edge[, 2] == tmp1[2])

which(test$edge[, 2] == 2)
which(test$edge[, 2] == 137648)

test$edge[c(1,2,137647), ]
sum(test$edge.length[c(1,2,137647)])

#distance is 0.00055

