## Function to identify and remove Species assignments for queries with close neighbors from different Species

## I think we can alter this function to work at any taxonomic level 
## and to work on the original reference tree tips...
## would be a better way to identify potentially ambigious tips than the current method I am using.


library(castor)
library(stringr)

nearest.neighbor.distances <- function(tax.df, placement.object, reference.tree, max.radius){
  
  #make sure all labels on the reference tree are unique
  reference.tree$label <- make.unique(reference.tree$label)
  
  #get the query sequence names
  queries.w.species <- tax.df %>% 
    filter(!is.na(Species)) %>% 
    pull(query.name)
  
  #generate return dataframe
  #which shows the query.name
  #the placement for that query
  #the light_weight_ratio
  #the nearest_dist to its neighbour
  dists <- data.frame("query.name"=NA, "Placement"=NA, "like_weight_ratio"=NA, "Nearest_Dist"=NA)
  #loop through each query sequence
  
  for(nm in queries.w.species){ 
    if(!"tbl_tree" %in% class(reference.tree)){class(reference.tree) <- c("tbl_tree", class(reference.tree))}
    plotTree <- as.phylo(reference.tree)
    ## Extract placements
    plc <- placement.object@placements[which(placement.object@placements$name==nm),] %>% arrange(desc(like_weight_ratio))
    plc <- plc %>% filter(like_weight_ratio > 0.5*max(like_weight_ratio))
    #we should be filtering plc by 0.5*max(like_weight_ratio)
    
    
    print(paste(which(queries.w.species==nm), "of", length(queries.w.species)))
    
    #go through each placement for that query sequence
    for(pind in 1:nrow(plc)){
      ## Add placement into tree
      plotTree <- as.phylo(reference.tree)
      plotTree <- bind.tip(plotTree, tip.label=paste0("Placement_", pind), 
                           edge.length=plc$pendant_length[pind], where=plc$node[pind], 
                           position=plc$distal_length[pind])
      
      plc1 <- which(plotTree$tip.label==paste0("Placement_", pind))
      
      ## Get neighbors within a given radius
      neighbors <- extract_tip_radius(plotTree, focal_tip = paste0("Placement_", pind), radius = max.radius, include_subtree = FALSE)
      
      ## if there are no other tips within the max.radius then we just return NA as the nearest_dist
      if(length(neighbors$radius_tips)==1){
        dists <- rbind(dists, c(nm, paste0("Placement_", pind), 
                                plc[pind,] %>% pull(like_weight_ratio), 
                                NA))
        break()
      } else {
        neighbors <- extract_tip_radius(plotTree, focal_tip = paste0("Placement_", pind), radius = max.radius, include_subtree = TRUE)
      }
      
      #grab taxa info from reference tree
      labelInfo <- reference.tree %>% select(label, Kingdom, Phylum, Class, Order, Family, Genus, Species)
      
      plotTree2 <- left_join(as_tibble(neighbors$subtree), labelInfo, by="label")
      
      plotTree2$isTip <- isTip(plotTree2, plotTree2$node)
      
      if(!"tbl_tree" %in% class(plotTree2)){class(plotTree2) <- c("tbl_tree", class(plotTree2))}
      
      #grab the current name of the placement before checking for differences in distances...
      parName <- tax.df %>% filter(query.name==nm) %>% pull(Species)
      
      #get the indexes that have different species labels than the par name
      diff.neighbors <-plotTree2  %>% 
        filter(is.na(Species) | !str_detect(string=parName, pattern=Species)) %>% 
        filter(isTip) %>% pull(node)
      
      
      placement_node <- plotTree2 %>% filter(label==paste0("Placement_", pind)) %>% pull(node)
      
      #we want to remove the node that the placement was at 
      #I think we can just do this by tracking the label rather than calling this function
      descendents_of_index <- offspring(neighbors$subtree, parent(neighbors$subtree, placement_node), type = "tips")
      
      
      #eligible.tips <- setdiff(diff.neighbors, placement_node) # change this as test
      eligible.tips <- setdiff(diff.neighbors, descendents_of_index)
      
      #If there are no tips within the max radius that have different species labels we return NA
      if(length(eligible.tips)==0){
        dists <- rbind(dists, c(nm, paste0("Placement_", pind), 
                                plc[pind,] %>% pull(like_weight_ratio), 
                                NA))} 
      #if there are tips left we find the nearest tip in the tree and grab its distance
      #save that distance to report back.
      else {      
            nearest_neighbor <- find_nearest_tips(neighbors$subtree, target_tips = eligible.tips)
            dists <- rbind(dists, c(nm, paste0("Placement_", pind), 
                                    plc[pind,] %>% pull(like_weight_ratio), 
                                    nearest_neighbor$nearest_distance_per_tip[placement_node]))
                                  
      }
    }
  }
  
  ## Format and calculate minimum distance per query
  dists <- dists %>% filter(!is.na(query.name))
  dists$like_weight_ratio <- as.numeric(dists$like_weight_ratio)
  dists$Nearest_Dist <- as.numeric(dists$Nearest_Dist)
  
  #why do we remove 2e-12 distances?
  
  #in this case we look at all the placements and just grab the near_dist as the minimum among all placements
  dists <- dists %>% 
    group_by(query.name) %>%
    filter(!is.na(Nearest_Dist)) %>%
    filter(Nearest_Dist!=2e-12) %>%
    summarize(minDist = min(Nearest_Dist, na.rm=T))
  
  return(dists)
  
}

nearest.neighbor.revisions <- function(tax.df, distances, radius){
  revise <- distances %>% filter(minDist<=radius ) %>% pull(query.name) 
  tax.revised <- tax.df %>%
    mutate(Species = ifelse(query.name %in% revise, NA, Species))
  result <- tax.revised
  return(result)
}
