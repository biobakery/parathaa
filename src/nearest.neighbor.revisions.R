## Function to identify and remove Species assignments for queries with close neighbors from different Species

library(castor)
library(stringr)

nearest.neighbor.distances <- function(tax.df, placement.object, reference.tree, max.radius){
  
  queries.w.species <- tax.df %>% 
    filter(!is.na(Species)) %>% 
    pull(query.name)
  
  dists <- data.frame("query.name"=NA, "Placement"=NA, "like_weight_ratio"=NA, "Nearest_Dist"=NA)
  for(nm in queries.w.species){ 
    if(!"tbl_tree" %in% class(reference.tree)){class(reference.tree) <- c("tbl_tree", class(reference.tree))}
    plotTree <- as.phylo(reference.tree)
    ## Extract placements
    plc <- placement.object@placements[which(placement.object@placements$name==nm),] %>% arrange(desc(like_weight_ratio))
    print(paste(which(queries.w.species==nm), "of", length(queries.w.species)))
    for(pind in 1:nrow(plc)){
      ## Add placement into tree
      plotTree <- as.phylo(reference.tree)
      plotTree <- bind.tip(plotTree, tip.label=paste0("Placement_", pind), 
                           edge.length=plc$pendant_length[pind], where=plc$node[pind], 
                           position=plc$distal_length[pind])
      plc1 <- which(plotTree$tip.label==paste0("Placement_", pind))
      
      ## Get neighbors within a given radius
      neighbors <- extract_tip_radius(plotTree, focal_tip = paste0("Placement_", pind), radius = max.radius, include_subtree = FALSE)
      if(length(neighbors$radius_tips)==1){
        dists <- rbind(dists, c(nm, paste0("Placement_", pind), 
                                plc[pind,] %>% pull(like_weight_ratio), 
                                NA))
        break()
      } else {
        neighbors <- extract_tip_radius(plotTree, focal_tip = paste0("Placement_", pind), radius = max.radius, include_subtree = TRUE)
      }
      
      labelInfo <- reference.tree %>% select(label, Kingdom, Phylum, Class, Order, Family, Genus, Species)
      plotTree2 <- left_join(as_tibble(neighbors$subtree), labelInfo, by="label")
      plotTree2$isTip <- isTip(plotTree2, plotTree2$node)
      if(!"tbl_tree" %in% class(plotTree2)){class(plotTree2) <- c("tbl_tree", class(plotTree2))}
      parName <- tax.df %>% filter(query.name==nm) %>% pull(Species)
      
      diff.neighbors <-plotTree2  %>% 
        filter(is.na(Species) | !str_detect(string=parName, pattern=Species)) %>% 
        filter(isTip) %>% pull(node)
      placement_node <- plotTree2 %>% filter(label==paste0("Placement_", pind)) %>% pull(node)
      descendents_of_index <- offspring(neighbors$subtree, parent(neighbors$subtree, placement_node), type = "tips")
      #eligible.tips <- setdiff(diff.neighbors, placement_node) # change this as test
      eligible.tips <- setdiff(diff.neighbors, descendents_of_index)
      
      if(length(eligible.tips)==0){
        dists <- rbind(dists, c(nm, paste0("Placement_", pind), 
                                plc[pind,] %>% pull(like_weight_ratio), 
                                NA))} else {      
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
