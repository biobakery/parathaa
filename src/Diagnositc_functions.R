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
library(stringr)

benchmark_classifications <- function(assignments, true_assignments, ret_frame=F, word_sep=" ", is_SILVA=TRUE){
  
  merged_data <- assignments %>% left_join(true_assignments, by="ID")
  
  if(is_SILVA){
    merged_data <- merged_data %>% 
      mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2, word_sep), collapse = ";" ))))
  }
  
  merged_data <- merged_data %>% 
    mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x),
           Genus.x = ifelse(Genus.x=="", NA, Genus.x))
  
  if(is_SILVA){
    merged_data2 <- merged_data %>% 
      dplyr::rowwise() %>%
      mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
             Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
      )
  }else{
    merged_data2 <- merged_data %>% 
      dplyr::rowwise() %>%
      mutate(Flag = ifelse(is.na(Species.x), NA, Species.y  %in% str_split(Species.x, ";", simplify = T)),
             Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
      )
  }
  
  
  
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
    return(list(results_specs, results_genus, merged_data2))
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


plot_placement <- function(thresholded_tree, level, placement, label, collapse=F, prop=0.01){
  ret_plots <- list()
  node_check <- paste("is",level,"Node", sep="")
  node_check <- sym(node_check)
  level <- sym(level)
  
  label_info <- placement@placements %>% filter(name==label)  %>% filter(like_weight_ratio==max(like_weight_ratio))
  
  thresholded_tree_phylo <- as.phylo(thresholded_tree)
  placed_node <- label_info$node

  for(i in placed_node){
    placed_at_tip=FALSE
    if(i <= length(thresholded_tree_phylo$tip.label)){
      #if given a tip we grab its parent node
      parent_node <- thresholded_tree$parent[thresholded_tree$node==i]
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
      placed_index <- which(keep_tips==i)
      sub_tips$PLACED_HERE[which(sub_tips$node==placed_index)] <- TRUE
    }
    
    keep_nodes <- c(subtree$new2old_node + length(thresholded_tree_phylo$tip.label))
    sub_nodes <- thresholded_tree[which(thresholded_tree$node %in% keep_nodes),]
    sub_nodes$node <- match(sub_nodes$node, (subtree$new2old_node  + length(thresholded_tree_phylo$tip.label))) +
      length(sub_tips$node)
    sub_nodes$PLACED_HERE <- FALSE
    if(!placed_at_tip){
      placed_index <- which(keep_nodes==i)
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
    ret_plots[[as.character(i)]] <- plot
    plot
  }

  return(ret_plots)
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
get_subtree_plot_data <- function(tree_df, node, isTip=FALSE, highlight_lab,
                                  level="Phylum"){
  
  tree_phylo <- as.phylo(tree_df)
  message("test")
  if(isTip){
    subtree <- get_subtree_at_node(tree_phylo, node)
    subtree_df <- as_tibble(subtree$subtree)

  }else{
    ntips <- length(which(tree_df$isTip))
    message("test")
    subtree <- get_subtree_at_node(tree_phylo, node=node-ntips)
    subtree_df <- as_tibble(subtree$subtree)
  }

  subtree_df$Phylum <- unlist(tree_df[,level])[match(subtree_df$label, tree_df$label)]
  subtree_df$highlight <- FALSE
  subtree_df$highlight[which(subtree_df$label %in% highlight_lab)] <- TRUE
  return(subtree_df)
}

get_n_parents <- function(tree_df, node, steps){
  i <- 0
  while(i < steps){
    i <- i + 1
    node <- tree_df$parent[which(tree_df$node==node)]
  }
  return(node)
}


#so for the jaccard distance with ambigiouty....
#whats the fastest way to calculate it...
#if something is x;y
#there is four possible scenarios
#either x is there
# y is there
# x and y is there (how do we deal with this one and do we care?)


# we can start with the simpleist implementation where if something is x;y we check if either x or y is in 
# if it is then we divide by the number of members and add?
# i'm not sure if that is the best way to do that but lets start with something easy

#if there is an ambigious that matches another ambigious
#we just divide the value being added by the total number of possible options...

## columns are samples and rows are taxon
## this calculation is going to be really slow how to speed it up...
#ambig_jaccard_distance <- function(matrix){
#   #first we need to remove rows that have ambigious assignments
#   #agreement
#   ret_matrix <- matrix(data=NA, nrow=ncol(matrix), ncol=ncol(matrix))
#   for(i in 1:ncol(matrix)){
#     for(j in 1:ncol(matrix)){
#       sample_i_taxon <- rownames(tmp_matrix)[which(tmp_matrix[,i] > 0)]
#       sample_j_taxon <- rownames(tmp_matrix)[which(tmp_matrix[,j]) > 0]
#       
#       #grab the non ambigious taxon
#       sample_i_taxon_non_ambig <- sample_i_taxon[which(grepl(";", sample_i_taxon, invert=TRUE))]
#       sample_j_taxon_non_ambig <- sample_j_taxon[which(grepl(";", sample_j_taxon, invert=TRUE))]
#       
#       #normal jaccard values without considering ambigious
#       reg_intersect <- intersect(sample_i_taxon_non_ambig, sample_j_taxon_non_ambig)
#       reg_union <- length(sample_i_taxon_non_ambig) + length(sample_j_taxon_non_ambig) - reg_intersect
#       
#       sample_i_total_taxon <- unlist(str_split(sample_i_taxon_non_ambig, ";"))
#       sample_j_total_taxon <- unlist(str_split(sample_j_taxon_non_ambig, ";"))
#       
#       ambig_taxon_i <- sample_i_taxon[which(grepl(";", sample_i_taxon, invert=FALSE))]
#       
#       ## so problem here becomes the double counting of j and i for amibigious taxon...
#       ## I think what we have to do is count the number of ambig shared
#       ## and then the number of unqieu ambig in each sample
#       ## add them all together and we get the dominator for the distance?
#       ## think about this more tomorrow....
#       ## if we then check the shared they are only counted once but the problem becomes what if an ambigious matches
#       ## an ambigious then we need to know the total number of taxon possible and divide hmmm
#       ## need to think about this more...
#       total_amib_value_i <- 0
#       total_amib_value_j <- 
#       for(ambig_taxon_i){
#         possible_taxon_i <- strsplit(ambig_taxon_i)
#         dom <- length(possible_taxon_i)
#         to_match <- unlist(possible_taxon_i)
#         if(grepl(paste(to_match, collapse="|"), sample_j_total_taxon)){
#           total_amib_value_i <- total_amib_value_i + 1/dom
#         }
#       }
#       for(amibig_taxon_j){
#         possible_taxon_j <- str_split(ambig_taxon_j)
#         dom <- length()
#       }
#       
#     }
#   }
#   return(ret_matrix)
# }

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection(a,b)
  return (intersection/union)
}
