#!/usr/bin/env Rscript
require(docopt)
'Usage:
   Exact_match.R [-q <query alignment file> -o <output directory> --threads <cores> -r <reference alighment file> -t <named tree file> --util1 <file containing species editor function>]

Options:
   -q query alignment file
   -o output directory
   --threads number of threads to run on
   -r reference alignment
   -t tree file
   --util1 script for species editting

 ]' -> doc
opts <- docopt(doc)

library(seqinr)
library(dplyr)
library(doSNOW)
parathaaDir <- (opts$p)
source(opts$util1)
taxa_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")


## function that will determine if there is a dominate taxa via binomial error model 
## if there are more than 1 exact match with differing taxonomies
binomial_check <- function(falseNegRate=0.05, acceptableProb=0.20, nodeGroups){
  
  if(pbinom(sum(nodeGroups)-max(nodeGroups)-1,
            sum(nodeGroups),falseNegRate, lower.tail = F) > acceptableProb){
    # if that is correct we then assign the int node to that phylum
    result <- names(nodeGroups)[which(nodeGroups==max(nodeGroups))][[1]]
  }else {result <- paste(
      names(which(nodeGroups >
                    qbinom(acceptableProb , sum(nodeGroups),falseNegRate, lower.tail = F))), collapse=";")
  }
  
}

#read in query alignment file
query_alignment <- read.alignment(opts$q, format="fasta")

#read in reference alignment file
ref_alignment <- read.alignment(opts$r, format="fasta")


## convert . into gaps
query_alignment[[3]] <- gsub("\\.", "-", query_alignment[[3]])
ref_alignment[[3]] <- gsub("\\.", "-", ref_alignment[[3]])

#fix ref names
ref_alignment[[2]] <- gsub("\t.*", "", ref_alignment[[2]])

#Get reference taxdata from SILVA:
load(opts$t)
tree <- resultData$tax_bestcuts

query_with_exact <- c()
#assignments <- data.frame(matrix(ncol=10, nrow=0))
#colnames(assignments) <- c("query.name", taxa_levels, "maxDist", "ref_matches")

cl <- makeCluster(as.numeric(opts$threads))
registerDoSNOW(cl)


assignments <- foreach(i=1:length(query_alignment[[2]]), .combine=bind_rows, .packages = c("dplyr", "seqinr")) %do% {
  query <- query_alignment[[3]][[i]]
  query_name <- query_alignment[[2]][[i]]
  
  #find exact matches between the ref_alighnment and the query
  matches <- which(unlist(ref_alignment[[3]]) %in% query)
  if(length(matches)!=0){
    ref_hit_id <- ref_alignment[[2]][matches]
    query_with_exact <- c(query_with_exact, query_name)
    
    ref_taxonomy <- tree %>% filter(label %in% ref_hit_id) %>%
      select(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
    
    #get number of unique assignments
    if(nrow(ref_taxonomy) > 1){
      ret_taxa <- c()
      for(level in taxa_levels){
        nodeGroups <- table(ref_taxonomy[,level])
        if(length(nodeGroups)==0){
          ret_taxa <- c(ret_taxa, NA)
        }else{
          ret_taxa <- c(ret_taxa, binomial_check(nodeGroups = nodeGroups))
        }

        
        temp_assignment <- data.frame(query.name=query_name,
                                      Kingdom=ret_taxa[1],
                                      Phylum=ret_taxa[2],
                                      Class=ret_taxa[3],
                                      Order=ret_taxa[4],
                                      Family=ret_taxa[5],
                                      Genus=ret_taxa[6],
                                      Species=ret_taxa[7],
                                      maxDist=0,
                                      ref_matches=paste0(ref_hit_id, collapse = ";"),
                                      row.names = "query.name")
      }
    }else{
      temp_assignment <- data.frame(query.name=query_name,
                                    Kingdom=ref_taxonomy$Kingdom,
                                    Phylum=ref_taxonomy$Phylum,
                                    Class=ref_taxonomy$Class,
                                    Order=ref_taxonomy$Order,
                                    Family=ref_taxonomy$Family,
                                    Genus=ref_taxonomy$Genus,
                                    Species=ref_taxonomy$Species,
                                    maxDist=0,
                                    ref_matches=paste0(ref_hit_id, collapse = ";"),
                                    row.names="query.name")
    }
    return(temp_assignment)
  }
}

## okay we need to remove lower level amibigious stuff if there is one in higher level.

for (i in 1:nrow(assignments)) {
  if(grepl(";", assignments[i,"Kingdom"])){
    assignments[i,taxa_levels[-1]] <- NA
  }else if(grepl(";", assignments[i, "Phylum"])){
    assignments[i, taxa_levels[-c(1,2)]] <- NA
  }else if(grepl(";", assignments[i, "Class"])){
    assignments[i, taxa_levels[-c(1,2,3)]] <- NA
  }else if(grepl(";", assignments[i, "Order"])){
    assignments[i, taxa_levels[-c(1,2,3,4)]] <- NA
  }else if(grepl(";", assignments[i, "Family"])){
    assignments[i, taxa_levels[-c(1,2,3,4,5)]] <- NA
  }else if(grepl(";", assignments[i, "Genus"])){
    assignments[i, taxa_levels[-c(1,2,3,4,5,6)]] <- NA
  }
}


## okay now we need to write out the new filtered alignment file

query_alignment_filt_names <- query_alignment[[2]][-which(query_alignment[[2]] %in% query_with_exact)]
query_alignment_filt_seqs <- query_alignment[[3]][match(query_alignment_filt_names, query_alignment[[2]])]

#write out the filtered alignment file that can be passed into parathaa's main algo.

write.fasta(sequences = query_alignment_filt_seqs, names=query_alignment_filt_names, file.out = paste0(opts$o, "/query_alignment_filt.fasta"))

## write assignments
write.table(assignments, file=file.path(opts$o, "taxonomic_assignments_exact.tsv"),sep = '\t', quote = F, row.names=F)
