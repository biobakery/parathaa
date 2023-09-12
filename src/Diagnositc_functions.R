### diagnosistic scripts...

#tax_file should be address to the taxonomic file
#region tree is address to newick file

#source("src/SILVA.species.editor.R")
library(dplyr)
library(forcats)
library(ggtree)
library(castor)
library(treeio)
library(tidyr)

benchmark_classifications <- function(assignments, true_assignments, ret_frame=F){
  
  merged_data <- assignments %>% left_join(true_assignments, by="ID")
  
  merged_data <- merged_data %>% 
    mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
  
  merged_data <- merged_data %>% 
    mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x),
           Genus.x = ifelse(Genus.x=="", NA, Genus.x))
  
  merged_data2 <- merged_data %>% 
    dplyr::rowwise() %>%
    mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
           Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
    )
  
  
  multi_spec <- merged_data2[grep(";", merged_data2$Species.x),]
  non_multi_spec <- merged_data2[grep(";", merged_data2$Species.x, invert = T), ]
  

  
  results_specs <- data.frame(
    unique_correct <- c(length(which(non_multi_spec$Flag)), 
                        length(which(non_multi_spec$Flag))/dim(assignments)[1]),
    one_to_many <- c(length(which(multi_spec$Flag)),
                     length(which(multi_spec$Flag))/dim(assignments)[1]),
    fpr <- c(length(which(!non_multi_spec$Flag)) + length(which(!multi_spec$Flag)),
             (length(which(!non_multi_spec$Flag)) + length(which(!multi_spec$Flag)))/dim(assignments)[1]),
    unassigned <- c(length(which(is.na(non_multi_spec$Flag))), 
                    length(which(is.na(non_multi_spec$Flag)))/dim(assignments)[1])
  )
  
  colnames(results_specs) <- c("unique_correct", "one_to_many", "fpr", "unassigned")

  multi_genera <- merged_data2[grep(";", merged_data2$Genus.x),]
  non_multi_genera <- merged_data2[grep(";", merged_data2$Genus.x, invert = T), ]
  
  results_genus <- data.frame(
    unique_correct <- c(length(which(non_multi_genera$Flag.genus)), 
                        length(which(non_multi_genera$Flag.genus))/dim(assignments)[1]),
    one_to_many <- c(length(which(multi_genera$Flag.genus)),
                     length(which(multi_genera$Flag.genus))/dim(assignments)[1]),
    fpr <- c(length(which(!non_multi_genera$Flag.genus)) + length(which(!multi_genera$Flag.genus)),
             (length(which(!non_multi_genera$Flag.genus)) + length(which(!multi_genera$Flag.genus)))/dim(assignments)[1]),
    unassigned <- c(length(which(is.na(non_multi_genera$Flag.genus))), 
                    length(which(is.na(non_multi_genera$Flag.genus)))/dim(assignments)[1])
  )
  
  colnames(results_genus) <- c("unique_correct", "one_to_many", "fpr", "unassigned")
  
  if(ret_frame)
    return(list(results_specs, results_genus))
  else
    return(list(results_specs, results_genus))
}

plot_region_specific_tree <- function(taxafile, region_tree, level="Phylum", prevelance=0.01,
                               layout="circular"){
  
  ## Bring in taxonomy file
  inFileTaxdata <- taxafile
  
  ## Tree made from database trimmed to region
  in.tree <- read.newick(region_tree)
  in.tree.data <- as_tibble(in.tree)
  
  if(grepl("\\|M", in.tree.data$label[1])){
    suppressWarnings({
      in.tree.data <- in.tree.data %>%
        separate(col=label, into=c("arbID", "primaryAccession"), remove=F, sep="\\|")
    })
  }else{
    suppressWarnings({
      in.tree.data <- in.tree.data %>%
        separate(col=label, into=c("primaryAccession", "arbID"), remove=F, sep="\\.")
    })
  }
  
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
  
  plot_data <- in.tree.data
  
  plot_data[,level] <- fct_lump(f=factor(unlist(plot_data[,level])), prop=prevelance)
  
  plot <- ggtree(in.tree, layout=layout) %<+% plot_data + geom_tippoint(aes(color=Phylum))
  
  return(plot)
}

plot_labeled_Int_nodes <- function(labeled_tree, level="Phylum", prevelance=0.01,
                                      layout="circular"){
  
  plot_data <- labeled_tree
  plot_data$has_level <- !is.na(plot_data[,level])
  labeled_tree_phylo <- as.phylo(labeled_tree)
  legend_label <- paste("Has", level, "Classification", sep=" ")
  plot <- ggtree(labeled_tree_phylo, layout=layout) %<+% plot_data + geom_nodepoint(aes(color=has_level), alpha=0.5)
  
  return(plot)
}

#given a taxonomic level and a tree
#it will calculate how well the agreement between that label and the tree is
#easiest way to do this is caluclate the average distance between members vs distance between a member and an outside member
#return that distance.
calculate_phylogenetic_taxonomic_agreement <- function(tree, level, taxdata, normalize=T){}

plot_placement <- function(thresholded_tree, level, placement, label, collapse=F, prop=0.01){
  node_check <- paste("is",level,"Node", sep="")
  node_check <- sym(node_check)
  level <- sym(level)
  
  placement <- get.placements(placement, by="best")
  label_info <- placement[which(placement$name==label),]  
  
  thresholded_tree_phylo <- as.phylo(thresholded_tree)
  placed_node <- label_info$node
  placed_at_tip=FALSE
  if(placed_node <= length(thresholded_tree_phylo$tip.label)){
    #if given a tip we grab its parent node
    parent_node <- thresholded_tree$parent[thresholded_tree$node==placed_node]
    subtree <- get_subtree_at_node(thresholded_tree_phylo, parent_node-length(thresholded_tree_phylo$tip.label))
    placed_at_tip=TRUE
  }else{
    subtree <- get_subtree_at_node(thresholded_tree_phylo, label_info$node-length(thresholded_tree_phylo$tip.label))
  }

  
  keep_tips <- subtree$new2old_tip
  sub_tips <- thresholded_tree[which(thresholded_tree$node %in% keep_tips),]
  sub_tips$node <- match(sub_tips$node, subtree$new2old_tip)
  sub_tips$PLACED_HERE <- FALSE
  if(placed_at_tip){
    placed_index <- which(keep_tips==placed_node)
    sub_tips$PLACED_HERE[which(sub_tips$node==placed_index)] <- TRUE
  }
  
  keep_nodes <- c(subtree$new2old_node + length(thresholded_tree_phylo$tip.label))
  sub_nodes <- thresholded_tree[which(thresholded_tree$node %in% keep_nodes),]
  sub_nodes$node <- match(sub_nodes$node, (subtree$new2old_node  + length(thresholded_tree_phylo$tip.label))) +
    length(sub_tips$node)
  sub_nodes$PLACED_HERE <- FALSE
  if(!placed_at_tip){
    placed_index <- which(keep_nodes==placed_node)
    sub_nodes$PLACED_HERE[which(sub_nodes$node==placed_index)] <- TRUE
  }
  
  thresholded_subtree <- rbind(sub_nodes, sub_tips)
  
  text_label1 <- paste("Distal length: ", round(label_info$distal_length, digits = 6))
  text_label2 <- paste("Pendant length: ", round(label_info$pendant_length, digits = 6))
  
  full_label <- paste(text_label1, text_label2)
  
  if(collapse){
    
    thresholded_subtree[,level] <- fct_lump(f=factor(unlist(thresholded_subtree[,level])), prop=prop)
  }
  
  if(placed_at_tip){
    plot <- ggtree(subtree$subtree) %<+% thresholded_subtree + geom_tippoint(aes(color=!!level, size=PLACED_HERE, 
                                                                                 shape=!!node_check)) +
      geom_nodepoint(aes(shape=!!node_check, size=!!node_check)) +
      scale_size_manual(values=c(1, 5)) + geom_treescale() +
      annotate(geom="text", x=-Inf, y=Inf, label=full_label, vjust=1, hjust=-.5)
    
  }else{
    plot <- ggtree(subtree$subtree) %<+% thresholded_subtree + geom_tippoint(aes(color=!!level, shape=!!node_check)) +
      geom_nodepoint(aes(shape=!!node_check, color=!!level)) +
      scale_size_manual(values=c(1, 5)) + geom_treescale()+
      annotate(geom="text", x=-Inf, y=Inf, label=full_label, vjust=1, hjust=-.5)
  }


  
  return(plot)
}

## function that takes in tree and then highlights the nodes that have 
## placements on them that are not resolved at X taxonomic level
## this will give an idea of whats going on there...

plot_low_resolved_placements <- function(taxafile, region_tree, level="Phylum", prevelance=0.01,
                                         layout="circular"){
  
  ## Bring in taxonomy file
  inFileTaxdata <- taxafile
  
  ## Tree made from database trimmed to region
  in.tree <- read.newick(region_tree)
  in.tree.data <- as_tibble(in.tree)
  
  if(grepl("\\|M", in.tree.data$label[1])){
    suppressWarnings({
      in.tree.data <- in.tree.data %>%
        separate(col=label, into=c("arbID", "primaryAccession"), remove=F, sep="\\|")
    })
  }else{
    suppressWarnings({
      in.tree.data <- in.tree.data %>%
        separate(col=label, into=c("primaryAccession", "arbID"), remove=F, sep="\\.")
    })
  }
  
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
  
  plot_data <- in.tree.data
  
  #load in placement data
  
  placments <- get.placements(read.jplace(jplace), by="best")
  
  
  
}


#given a node plots the subtree
get_subtree_plot_data <- function(tree_df, node){
  
  tree_phylo <- as.phylo(tree_df)
  
  subtree <- get_subtree_at_node(tree_phylo, node)
  subtree_df <- as_tibble(subtree$subtree)
  
  subtree_df$Phylum <- tree_df$Phylum[match(subtree_df$label, tree_df$label)]
  
  return(subtree_df)
}
