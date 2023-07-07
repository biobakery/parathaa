## Compare taxonomies from DADA2 vs PARATHAA

library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(dplyr) # data handling 
library(stringr)
library(ggtree)
library(treeio)
library(dada2)
suppressPackageStartupMessages(library(seqinr))
#source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

##################################
## Synthetic community analysis ##
##################################

## Specify input paths/working directory
setwd("C:/Users/mshort/Documents/proj/parathaa/")
parathaaFile <- "C:/Users/mshort/Documents/proj/parathaa/output/20230526_SynthV1V2_weights3to1/taxonomic_assignments.tsv"
sequenceFile <- "input/SILVAsubsample_SeedGenera_V1V2.pcr.fasta"
DADAdb <- "input/20230111.silva.seed_v138_1.ng.dada.fasta"
DADAdb.sp <- "input/silva_species_assignment_v138.1.seedonly.fa"

# Read in parathaa data
parathaaData <- read.delim(
  parathaaFile,
  sep='\t', fill=T, stringsAsFactors = F, header=T)
tax_parathaa <- parathaaData %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  group_by(query.name) 

taxmat <- tax_parathaa %>% as.matrix
rownames(taxmat) <- tax_parathaa$query.name
taxmat <- taxmat[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa <- tax_table(taxmat)

otutab <- tax_parathaa %>% 
  dplyr::group_by(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(query.name, name="Parathaa") 
otumat <- as.matrix(otutab$Parathaa)
rownames(otumat) <- otutab$query.name
colnames(otumat) <- "Parathaa.V1V2"
OTU_parathaa <- otu_table(otumat, taxa_are_rows = TRUE)

samp_parathaa <- data.frame(colnames(OTU_parathaa), rep("parathaa", length(colnames(OTU_parathaa))),
                            rep("V1V2", length(colnames(OTU_parathaa))))
colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_parathaa) <- samp_parathaa$sampleID
str(samp_parathaa)
SAMP_parathaa <- sample_data(samp_parathaa)

ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)
print(ps1_parathaa)

## Assign taxonomy with DADA2
## First, get names and sequences from fasta file
getNames <- read.fasta(file(sequenceFile), as.string = TRUE,
                       forceDNAtolower = FALSE, whole.header = FALSE)
names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
names1 <- names1[,1] %>%
  str_remove(">")
name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)

## Next, assign taxonomy to genus level with DADA2 (takes a few minutes)
set.seed(3874)
taxa <- assignTaxonomy(sequenceFile, 
                       DADAdb,
                       multithread=TRUE)

## Remove sequences with undefined ("N") bases, store until after species assignment
taxa.test <- as.data.frame(taxa)
taxa.test$taxaIDs <- names1
nChars <- grep("N", rownames(taxa.test))
print(paste("Removing", length(nChars), "sequences with N bases"))
withNbases <- taxa.test[nChars,]
taxa <- taxa[-nChars,]

## Perform species assignment with DADA2 (takes a few minutes)
taxa.sp <- addSpecies(taxa, DADAdb.sp)

## Add in reference IDs and taxonomy from sequences with "N" bases
tax_dada <- as.data.frame(taxa.sp)
tax_dada$sequence <- str_split(rownames(tax_dada), "\\.", simplify=TRUE)[,1]
getnamSubset <-name.df %>% filter(sequence %in% tax_dada$sequence)
tax_dada2 <- cbind(tax_dada, "taxaIDs" =getnamSubset$taxaIDs)
tax_dada2 <- full_join(tax_dada2, withNbases)
rownames(tax_dada2) <- tax_dada2[,"taxaIDs"]

tax_dada3 <- tax_dada2 %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  as.matrix()
## Rename species to include genus names
tax_dada3 <- cbind(tax_dada3, "Species2"=NA)
tax_dada3[which(!is.na(tax_dada3[,"Species"])), "Species2"] <-  
  paste(tax_dada3[which(!is.na(tax_dada3[,"Species"])), "Genus"], 
        tax_dada3[which(!is.na(tax_dada3[,"Species"])), "Species"] )
tax_dada3 <- tax_dada3[,!colnames(tax_dada3) %in% "Species"]
colnames(tax_dada3)[which(colnames(tax_dada3)=="Species2")] <- "Species"

## Place DADA2 taxonomies into phyloseq object
TAX_dada <- tax_table(tax_dada3)

otutab <- as.data.frame(tax_dada2) %>% 
  select(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::group_by(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(taxaIDs, name="DADA2") 
otumat <- as.matrix(otutab$DADA2)
rownames(otumat) <- otutab$taxaIDs
colnames(otumat) <- "DADA2.V1V2"
OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)

samp_dada <- data.frame(colnames(OTU_dada), rep("DADA2", length(colnames(OTU_dada))),
                        rep("V1V2", length(colnames(OTU_dada))))
colnames(samp_dada) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_dada) <- samp_dada$sampleID
str(samp_dada)
SAMP_dada <- sample_data(samp_dada)

ps1_dada <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)

print(ps1_dada)


############################################
### Assess performance on Synthetic Data ###
############################################

#Get taxdata from Synthetic.Reads.R
source("src/Synthetic.Reads.R")

synth.parathaa<- as.data.frame(tax_table(ps1_parathaa))
synth.parathaa$AccID <- rownames(synth.parathaa)

synth.parathaa2 <- left_join(synth.parathaa, taxdata, by="AccID")
synth.parathaa2 <- synth.parathaa2 %>% 
  mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
synth.parathaa2 <- synth.parathaa2 %>% 
  mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x),
         Genus.x = ifelse(Genus.x=="", NA, Genus.x))
synth.parathaa3 <- synth.parathaa2 %>% 
  dplyr::rowwise() %>%
  mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
         Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
  )

synth.dada <- as.data.frame(tax_table(ps1_dada))
synth.dada$AccID <- rownames(synth.dada)

synth.dada2 <- left_join(synth.dada, taxdata, by="AccID")
synth.dada2 <- synth.dada2 %>% mutate(Species.x = word(Species.x, 1, 2))
synth.dada3 <- synth.dada2 %>% 
  rowwise() %>%
  mutate(Flag = ifelse(is.na(Species.x), 
                       NA, 
                       word(Species.y, 1, 2) %in% t(apply(str_split(Species.x, ";", simplify=T), 1, FUN=function(X) word(X, 1,2)))),
         Flag.genus = ifelse(is.na(Genus.x), 
                             NA, 
                             Genus.y %in% (str_split(Genus.x, ";", simplify=T))
         )
  )

compare.synth <- dplyr::full_join(synth.dada3, synth.parathaa3, by="AccID")
compare.synth <- compare.synth %>% 
  mutate(Species.silva = word(Species.y.y, 1, 2)) %>%
  rename(Species.parathaa = Species.x.y,
         Species.dada = Species.x.x,
         Genus.dada = Genus.x.x,
         Genus.parathaa = Genus.x.y,
         Genus.silva = Genus.y.y)


genus.summary <- data.frame("mWeight"=3, "sWeight"=1, "binErr"=0.05, "binP"=0.20)
genus.summary$unique.correct <- round(nrow(compare.synth %>% filter(Genus.parathaa==Genus.silva))/nrow(compare.synth), 3)
genus.summary$one.to.many <- round(nrow(compare.synth %>% filter(Flag.genus.y & Genus.parathaa!=Genus.silva))/nrow(compare.synth), 3)
genus.summary$fpr <- round(nrow(compare.synth %>% filter(!Flag.genus.y ))/nrow(compare.synth), 3)
genus.summary$unassigned <- round(nrow(compare.synth %>% filter(is.na(Flag.genus.y) ))/nrow(compare.synth), 3)
print(genus.summary)

species.summary <- data.frame("mWeight"=3, "sWeight"=1, "binErr"=0.05, "binP"=0.20)
species.summary$unique.correct <- round(nrow(compare.synth %>% filter(Species.parathaa==Species.silva))/nrow(compare.synth), 3)
species.summary$one.to.many <- round(nrow(compare.synth %>% filter(Flag.y & Species.parathaa!=Species.silva))/nrow(compare.synth), 3)
species.summary$fpr <- round(nrow(compare.synth %>% filter(!Flag.y ))/nrow(compare.synth), 3)
species.summary$unassigned <- round(nrow(compare.synth %>% filter(is.na(Flag.y) ))/nrow(compare.synth), 3)
print(species.summary)



compare.genus <- data.frame()
## Both uniquely correct:
compare.genus["Dada2_correct", "Parathaa_correct"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
## Paratha uniquely correct, dada2 incorrect:
compare.genus["Dada2_incorrect", "Parathaa_correct"] <- nrow(compare.synth %>% filter(!Flag.genus.x & !is.na(Genus.dada)  & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
## Paratha uniquely correct, dada2 unassigned:
compare.genus["Dada2_unassigned", "Parathaa_correct"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x)  & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

##Paratha partly correct, dada2 uniquely correct:
compare.genus["Dada2_correct", "Parathaa_1toMany"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha partly correct, dada2 incorrect:
compare.genus["Dada2_incorrect", "Parathaa_1toMany"] <- nrow(compare.synth %>% filter(!Flag.genus.x & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha partly correct, dada2 unassigned:
compare.genus["Dada2_unassigned", "Parathaa_1toMany"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x) & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

##Paratha incorrect, dada2 uniquely correct:
compare.genus["Dada2_correct", "Parathaa_incorrect"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Both incorrect:
compare.genus["Dada2_incorrect", "Parathaa_incorrect"] <- nrow(compare.synth %>% filter(!Flag.genus.x & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha incorrect, dada2 unassigned:
compare.genus["Dada2_unassigned", "Parathaa_incorrect"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x) & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

##Paratha unassigned, dada2 uniquely correct:
compare.genus["Dada2_correct", "Parathaa_unassigned"] <- nrow(compare.synth %>% filter(Genus.dada==Genus.silva & is.na(Flag.genus.y) ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha unassigned, dada2 incorrect:
compare.genus["Dada2_incorrect", "Parathaa_unassigned"] <- nrow(compare.synth %>% filter(!Flag.genus.x & is.na(Flag.genus.y)) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Both unassigned:
compare.genus["Dada2_unassigned", "Parathaa_unassigned"] <- nrow(compare.synth %>% filter(is.na(Flag.genus.x) & is.na(Flag.genus.y) ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

sum(compare.genus)
compare.genus$Sum <- apply(compare.genus, 1, sum)
compare.genus["Sum",] <- apply(compare.genus, 2, sum)

compare.genus[,"Pct"] <- round(compare.genus[,"Sum"] / compare.genus["Sum", "Sum"], 3)
compare.genus["Pct",] <- round(compare.genus["Sum",] / compare.genus["Sum", "Sum"], 3)

print("Genus comparison:")
print(t(compare.genus))

###################################
#### Species-level Summary ########
###################################

## Distinguish between query species "included" in the database and those "novel" to the database
taxdata_seed <- taxdata %>% 
  filter(primaryAccession %in% SeedTax$primaryAccession)
taxdata_SP <- word(taxdata_seed$Species, 1, 2)

compare.synth.full <- compare.synth
compare.synth.novel <- compare.synth.full %>% filter(!Species.silva %in% taxdata_SP)
compare.synth.included <- compare.synth.full %>% filter(Species.silva %in% taxdata_SP)

sp.list <- list(compare.synth.full, compare.synth.novel, compare.synth.included)
names(sp.list) <- c("All queries", "Novel Queries Only", "Included Queries Only")

for(ind in 1:length(sp.list)){
  comp1 <- sp.list[[ind]]
  compare.species <- data.frame()
  ## Both uniquely correct:
  compare.species["Dada2_correct", "Parathaa_correct"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ## Paratha uniquely correct, dada2 incorrect:
  compare.species["Dada2_incorrect", "Parathaa_correct"] <- nrow(comp1 %>% filter(!Flag.x & !is.na(Species.dada)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ## Paratha uniquely correct, dada2 unassigned:
  compare.species["Dada2_unassigned", "Parathaa_correct"] <- nrow(comp1 %>% filter(is.na(Flag.x)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  
  ##Paratha partly correct, dada2 uniquely correct:
  compare.species["Dada2_correct", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ##Paratha partly correct, dada2 incorrect:
  compare.species["Dada2_incorrect", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(!Flag.x & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ##Paratha partly correct, dada2 unassigned:
  compare.species["Dada2_unassigned", "Parathaa_1toMany"] <- nrow(comp1 %>% filter(is.na(Flag.x) & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  
  ##Paratha incorrect, dada2 uniquely correct:
  compare.species["Dada2_correct", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ##Both incorrect:
  compare.species["Dada2_incorrect", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(!Flag.x & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ##Paratha incorrect, dada2 unassigned:
  compare.species["Dada2_unassigned", "Parathaa_incorrect"] <- nrow(comp1 %>% filter(is.na(Flag.x) & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  
  ##Paratha unassigned, dada2 uniquely correct:
  compare.species["Dada2_correct", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(Species.dada==Species.silva & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ##Paratha unassigned, dada2 incorrect:
  compare.species["Dada2_incorrect", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(!Flag.x & is.na(Flag.y)) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  ##Both unassigned:
  compare.species["Dada2_unassigned", "Parathaa_unassigned"] <- nrow(comp1 %>% filter(is.na(Flag.x) & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
  
  sum(compare.species)
  compare.species$Sum <- apply(compare.species, 1, sum)
  compare.species["Sum",] <- apply(compare.species, 2, sum)
  
  compare.species[,"Pct"] <- round(compare.species[,"Sum"] / compare.species["Sum", "Sum"], 3)
  compare.species["Pct",] <- round(compare.species["Sum",] / compare.species["Sum", "Sum"], 3)
  
  print(names(sp.list)[ind])
  print(t(compare.species))
  
}


###################################
##    Examine Relevant Trees     ##
## Not a clean script but useful ##
###################################

#save(compare.synth.full, file = "output/20230120_SyntheticV1V2_removeUndefSp/compare.synth.RData")

## Load in Data
load("output/20230120_SyntheticV1V2_removeUndefSp/resultTree_bestThresholds.RData")
load("output/20230120_SyntheticV1V2_removeUndefSp/compare.synth.RData")
in.jplace <- read.jplace("output/20230120_SyntheticV1V2_removeUndefSp/merged.sub.jplace")


in.tree.data <- resultData$tax_bestcuts
## Example: "KT962913.1.1234" is a short sequence (27 NT) that dada2 assigns correctly to genus level but parathaa has wrong phylum
# Look at taxonomic assignments for a given query
compare.synth.full %>% filter(AccID=="AB930131.1.1469") %>% View
plotTree <- tree_subset(as.treedata(in.tree.data), node = 7, levels_back = 10)
ggtree(plotTree, aes(color=Genus)) + geom_tippoint() + geom_nodepoint()

## Identify queries based on how they're assigned
#Dada2 right, parathaa wrong
ids <- compare.synth %>% filter(Species.dada==Species.silva & !Flag.y ) %>% select(AccID) %>% as.data.frame()
ids <- ids$AccID
#Dada2 unassigned, paratha wrong
ids <- compare.synth %>% filter(is.na(Flag.x) & !Flag.y ) %>% select(AccID) %>% as.data.frame()
ids <- ids$AccID
# Both correct 
ids <- compare.synth.full %>% filter(Species.dada==Species.silva & Species.parathaa==Species.silva)
ids <- ids$AccID
## Dada2 Genus right, parathaa wrong:
ids <- compare.synth %>% filter(Genus.dada==Genus.silva & !Flag.genus.y ) %>% select(AccID) %>% as.data.frame()
ids <- ids$AccID

#Choose which set of ids to look at from above, plot trees for a few instances
for(nm in ids[1:20]){ # 
  plotTree <- tree_subset(as.treedata(in.tree.data), node = in.jplace@placements$node[which(in.jplace@placements$name==nm)][1], levels_back = 2)
  truName <- compare.synth %>% filter(AccID==nm) %>% select(Genus.silva) %>% as.character()
  parName <- compare.synth %>% filter(AccID==nm) %>% select(Genus.parathaa) %>% as.character()
  ggtree(plotTree, aes(color=Genus)) + geom_tippoint() + geom_nodepoint() + geom_tiplab(aes(label=AccID)) + labs(title=paste0("Source: ", truName), 
                                                                                                                 subtitle = paste0("Parathaa: ", parName))
  ggsave(paste0("output/20230120_SyntheticV1V2_removeUndefSp/genusDadacorrect/", nm, ".png"))
  
}

#Look at placements for a given query
in.jplace@placements %>% filter(name=="AB009013.1.1439") %>% View
