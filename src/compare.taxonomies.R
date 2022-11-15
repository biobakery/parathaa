## Compare taxonomies from SILVA vs PARATHAA

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(tidyverse) # data handling  
library(gridExtra)
library(lemon)
library(vegan)
library(microbiomeutilities)
library(pheatmap)
library(dada2)
library (devtools)
source_url("https://raw.githubusercontent.com/lrjoshi/FastaTabular/master/fasta_and_tabular.R")

## Read in data with taxonomy from Kelsey's analysis with SILVA
kelseyDataRaw <- read.delim("~/proj/PPITAA/input/all_samples_taxonomy_closed_reference.tsv", 
                             sep='\t', fill=T, stringsAsFactors = F, header=T)
tax_oldDADA2 <- kelseyDataRaw %>%
  dplyr::select(taxonomy) %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family",
                            "Genus", "Species"), sep=";")

## Remove samples as instructed by Kelsey
kelseyDataFull <- kelseyDataRaw %>%
  dplyr::select(-c(X,
            CBSF036_S66_L001, 
            CBSF060_S82_L001, 
            CBSF078_S77_L001, 
            CBSF079_S65_L001,
            Undetermined_S0_L001,
            taxonomy))

## Full SILVA DADA2 training set
#taxa <- assignTaxonomy("~/proj/PPITAA/input/amplicons_fromKelsey.fasta", 
#                       "/Users/mis696/proj/PPITAA/input/silva_nr99_v138.1_train_set.fa.gz", 
#                       multithread=TRUE)

## SILVA Seed-based DADA2 training set
taxa <- assignTaxonomy("~/proj/PPITAA/input/amplicons_fromKelsey.fasta", 
                       "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta", 
                       multithread=TRUE)
taxa.sp <- addSpecies(taxa, "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")


tax_dada <- as.matrix(taxa.sp)
#Get IDs from format_amplicon_table.R "kelseyData" object
rownames(tax_dada) <- paste0(kelseyData$taxaIDs, "a")

TAX_dada <- tax_table(tax_dada)

otumat <- as.matrix(kelseyDataFull)
rownames(otumat) <- paste0(kelseyData$taxaIDs, "a")
colnames(otumat) <- paste0(colnames(otumat), "a")
OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)

samp_dada <- data.frame(colnames(OTU_dada), rep("DADA2", length(colnames(OTU_dada))))

colnames(samp_dada) <- c("sampleID", "Taxonomy_type")
rownames(samp_dada) <- samp_dada$sampleID
str(samp_dada)
SAMP_dada <- sample_data(samp_dada)

ps1_dada <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)

print(ps1_dada)



## Read in PARATHAA results
parathaData <- read.delim("/Users/mis696/proj/PPITAA/output/20220926_v4v5/taxonomic_assignments.tsv", 
                          sep='\t', fill=T, stringsAsFactors = F, header=T)
parathaData$query.num <- as.numeric(gsub("ID", "", parathaData$query.name))
parathaData <- parathaData %>% arrange(query.num)

tax_paratha <- parathaData %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  group_by(query.name) %>% 
  filter(row_number()==1)
#tax$Kingdom[which(!is.na(tax$Kingdom))] <- paste0("k__", tax$Kingdom[which(!is.na(tax$Kingdom))])
#tax$Phylum[which(!is.na(tax$Phylum))] <- paste0("p__", tax$Phylum[which(!is.na(tax$Phylum))])
#tax$Class[which(!is.na(tax$Class))] <- paste0("c__", tax$Class[which(!is.na(tax$Class))])
#tax$Order[which(!is.na(tax$Order))] <- paste0("o__", tax$Order[which(!is.na(tax$Order))])
#tax$Family[which(!is.na(tax$Family))] <- paste0("f__", tax$Family[which(!is.na(tax$Family))])
#tax$Genus[which(!is.na(tax$Genus))] <- paste0("g__", tax$Genus[which(!is.na(tax$Genus))])
#tax$Species[which(!is.na(tax$Species))] <- paste0("s__", tax$Species[which(!is.na(tax$Species))])

taxmat_paratha <- as.matrix(tax_paratha)
rownames(taxmat_paratha) <- paste0(tax_paratha$query.name, "b")
taxmat_paratha <- taxmat_paratha[,-1]
TAX_parathaa <- tax_table(taxmat_paratha)

otumat_parathaa <- kelseyDataFull
otumat_parathaa <- as.matrix(otumat_parathaa)
samp_parathaa <- data.frame(colnames(otumat_parathaa), rep("parathaa", length(colnames(otumat_parathaa))))
rownames(otumat_parathaa) <- paste0(kelseyData$taxaIDs, "b")
colnames(otumat_parathaa) <- paste0(colnames(otumat_parathaa), "b")
OTU_parathaa <- otu_table(otumat_parathaa, taxa_are_rows = TRUE)


colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type")
rownames(samp_parathaa) <- paste0(samp_parathaa$sampleID, "b")
str(samp_parathaa)
SAMP_parathaa <- sample_data(samp_parathaa)

ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)


ps1_all<- merge_phyloseq(ps1_dada, ps1_parathaa)
ps1.com <- ps1_all
#ps1.com <- ps1_parathaa

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

#getPalette = colorRampPalette(brewer.pal(colourCount, ))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object

# now edit the unclassified taxa
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Phylum"]), "Phylum"] <- "p__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Class"]), "Class"] <- "c__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Order"]), "Order"] <- "o__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Family"]), "Family"] <- "f__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Genus"]), "Genus"] <- "g__"
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Species"]), "Species"] <- "s__"
tax_table(ps1.com)[which(tax_table(ps1.com)[, "Species"]=="uncultured bacterium"), "Species"] <- "s__"

# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Now we need to plot at family level, we can do it as follows:

# first remove the phy_tree

ps1.com@phy_tree <- NULL


Levels <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

# Taxonomy plots
for(level in Levels){
  ps1.com.lev <- aggregate_rare(ps1.com, level, detection = .01/100, prevalence = 10/100)
  #ps1.com.lev <- ps1.com
  ps1.com.lev.agg <- aggregate_top_taxa2(ps1.com.lev, 10, level)
  #ps1.com.lev.agg <- aggregate_taxa(ps1.com.lev, level)
  ps1.com.lev.agg.rel <- microbiome::transform(ps1.com.lev.agg, "compositional")
  plot.composition.relAbun <- plot_composition(ps1.com.lev.agg.rel,
                                               #sample.sort = "sampleID",
                                               sample.sort = "Taxonomy_type",
                                               #otu.sort = "abundance",
                                               x.label = "SampleID") 
  xlabs <- rep("", nrow(sample_data(ps1.com.lev.agg.rel)))
  xlabs[ceiling(nrow(sample_data(ps1.com.lev.agg.rel))/4)] <- "DADA2"
  xlabs[ceiling(nrow(sample_data(ps1.com.lev.agg.rel))*3/4)+1] <- "parathaa"
  plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
  plot.composition.relAbun <- plot.composition.relAbun +  theme_bw()  + scale_fill_brewer(palette="Paired")  #scale_fill_manual(level, values = glasbey()) #
  #plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
  plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + #guide_italics + theme(legend.title = element_text(size = 18)) + 
    guides(fill= 
             guide_legend(
               ncol=1,
               label.theme = element_text(face="italic", size=14),
               title.theme = element_text(size=18),
               )
           )
  plot.composition.relAbun <- plot.composition.relAbun +  scale_x_discrete(labels=xlabs) + theme(axis.text.x = element_text(size=14, angle=0))
  #ggsave(filename=paste0("~/proj/parathaa/output/viz/", level, "parathaa_taxonomy.png"), plot.composition.relAbun, width=10, height=6)
  if(level!="Species")
    ggsave(filename=paste0("/Users/mis696/proj/parathaa/output/20221115_ASD_dadaseed", level, "taxonomy.png"), plot.composition.relAbun, width=10, height=6)
  if(level=="Species")
    ggsave(filename=paste0("/Users/mis696/proj/parathaa/output/20221115_ASD_dadaseed", level, "taxonomy.png"), plot.composition.relAbun, width=14, height=6)
}


### Ordination plots
plotList <- list()
for(level in Levels){
  ps1.com.lev <- aggregate_rare(ps1.com, level, detection = .01/100, prevalence = 10/100)
  ps1.com.lev.agg <- aggregate_taxa(ps1.com.lev, level=level)
  ps1.com.lev.agg.rel <- microbiome::transform(ps1.com.lev.agg, "compositional")
  GP.ord <- ordinate(ps1.com.lev.agg.rel, "PCoA", "bray")
  #brays <- phyloseq::distance(ps1.com.lev.agg.rel, method="bray")
  #adonis2(brays ~ sample_data(ps1.com.lev.agg.rel)$Taxonomy_type)
  p1 = plot_ordination(ps1.com.lev.agg.rel, GP.ord, type="samples", color="Taxonomy_type", title=level) + 
    geom_polygon(aes(group=sampleID), color="grey")
  plotList[[level]] <- p1
}
legend <- g_legend(plotList[["Phylum"]] + theme(legend.position='bottom'))
p1 <- grid.arrange(plotList[["Phylum"]]+theme(legend.position='hidden'), plotList[["Class"]]+theme(legend.position='hidden'), 
             plotList[["Order"]]+theme(legend.position='hidden'), plotList[["Family"]]+theme(legend.position='hidden'),
             plotList[["Genus"]]+theme(legend.position='hidden'), legend, nrow = 3)

ggsave(filename=paste0("~/proj/parathaa/output/viz/AllLevels_PCoA_bray_seedDada.png"), plot=p1, height = 8, width=6, dpi=300)



#Heatmap
dat1 <- bind_cols(data.frame(tax_table(ps1_parathaa))$Genus, data.frame(tax_table(ps1_dada))$Genus)
colnames(dat1) <- c("paratha_Genus", "dada2_Genus")

topBugsDada2 <- dat1 %>% group_by(dada2_Genus) %>% dplyr::count() 
topBugsDada2 <- topBugsDada2 %>% arrange(desc(n)) 
topBugsDada2Names <- data.frame(topBugsDada2[1:30,])$dada2_Genus


dat2 <- dat1 %>% 
  group_by(paratha_Genus, dada2_Genus) %>%
  replace_na(list(paratha_Genus = "Unknown", dada2_Genus = "Unknown")) %>%
  dplyr::count() %>% 
  filter(dada2_Genus %in% topBugsDada2Names) 

dat3 <- dat2 %>%
  pivot_wider(names_from = paratha_Genus, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "dada2_Genus")

pheatmap(dat3)#, scale="row")


###############################
### Mock community analysis ###
###############################

## read in fasta
taxa <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SRR3225703.fasta", 
                       "/Users/mis696/proj/PPITAA/input/silva_nr99_v138.1_train_set.fa.gz",
                       multithread=TRUE)
nChars <- grep("N", rownames(taxa))
print(paste("Removing", length(nChars), "sequences with N bases"))
taxa <- taxa[-nChars,]
taxa.sp <- addSpecies(taxa, "/Users/mis696/proj/PPITAA/input/silva_species_assignment_v138.1.fa.gz")
taxa.sp <- taxa

##taxa.sp2 <- taxa
otutab <- cbind(rownames(taxa.sp), taxa.sp)
otutab <- as_tibble(otutab) %>% dplyr::rename("ASV"="V1") %>% dplyr::count(ASV, name="SRR3225703") 
otumat <- as.matrix(otutab$SRR3225703)
rownames(otumat) <- otutab$ASV
colnames(otumat) <- "SRR3225703"
OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)


taxtab <- cbind(rownames(taxa.sp), taxa.sp)
taxtab <- as_tibble(taxtab) %>% 
  dplyr::rename("ASV"="V1") %>% 
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225703") 
taxmat <- taxtab %>% as.matrix
rownames(taxmat) <- taxtab$ASV
taxmat <- taxmat[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
taxmat <- cbind(taxmat, "Species2"=NA)
taxmat[which(!is.na(taxmat[,"Species"])), "Species2"] <-  
  paste(taxmat[which(!is.na(taxmat[,"Species"])), "Genus"], 
        taxmat[which(!is.na(taxmat[,"Species"])), "Species"] )
taxmat <- taxmat[,!colnames(taxmat) %in% "Species"]
colnames(taxmat)[which(colnames(taxmat)=="Species2")] <- "Species"

TAX_dada <- tax_table(taxmat)
#One sequence has 2 possible assignments, need to reassign taxon name
#rownames(TAX_dada)[16506] <- paste0(rownames(TAX_dada)[16506], ".1")


samp_dada <- data.frame(colnames(OTU_dada), rep("DADA2", length(colnames(OTU_dada))))
colnames(samp_dada) <- c("sampleID", "Taxonomy_type")
rownames(samp_dada) <- samp_dada$sampleID
str(samp_dada)
SAMP_dada <- sample_data(samp_dada)

ps1_dada <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)

## Read in parathaa data
parathaData <- read.delim("/Users/mis696/proj/parathaa/output/20221101_MiSeqV4V5Mock/taxonomic_assignments.tsv", 
                          sep='\t', fill=T, stringsAsFactors = F, header=T)
tax_paratha <- parathaData %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  group_by(query.name) %>% 
  filter(row_number()==1)

suppressPackageStartupMessages(library(seqinr))
parsed = seqinr::read.fasta(file('/Users/mis696/proj/parathaa/input/SRR3225703.fasta'), as.string = TRUE,
                    forceDNAtolower = FALSE, whole.header = FALSE)
table = data.frame("ASV"=unlist(parsed), "query.name" = sapply(parsed, attr, 'name'), row.names=NULL)
#View(tax_paratha)

tax_paratha <- merge(tax_paratha, table)  %>%
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225703") 

taxmat <- tax_paratha %>% as.matrix
rownames(taxmat) <- paste0("x", tax_paratha$ASV )
taxmat <- taxmat[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa <- tax_table(taxmat)


otumat <- as.matrix(tax_paratha$SRR3225703)
rownames(otumat) <- paste0("x", tax_paratha$ASV )
colnames(otumat) <- "SRR3225703b"
OTU_parathaa <- otu_table(otumat, taxa_are_rows = TRUE)

samp_parathaa <- data.frame(colnames(OTU_parathaa), rep("parathaa", length(colnames(OTU_parathaa))))
colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type")
rownames(samp_parathaa) <- samp_parathaa$sampleID
str(samp_parathaa)
SAMP_parathaa <- sample_data(samp_parathaa)

ps1_parathaa <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)

##Trying to make a dummy second sample to deal with aggregate_taxa issue

ps1_all<- merge_phyloseq(ps1_dada, ps1_parathaa)




## Heatmap ##
taxtab <- cbind(rownames(taxa.sp), taxa.sp)
taxtab <- as_tibble(taxtab) %>% 
  dplyr::rename("ASV"="V1") 

taxmat <- cbind(taxtab, "Species2"=NA)
taxmat[which(!is.na(taxmat[,"Species"])), "Species2"] <-  
  paste(taxmat[which(!is.na(taxmat[,"Species"])), "Genus"], 
        taxmat[which(!is.na(taxmat[,"Species"])), "Species"] )
taxmat <- taxmat[,!colnames(taxmat) %in% "Species"]
colnames(taxmat)[which(colnames(taxmat)=="Species2")] <- "Species"

taxmat <- taxmat %>%
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225703") 
    

dat1 <- merge(tax_paratha, taxmat, by="ASV")
dat1 <- dat1 %>%
  dplyr::rename("paratha_Genus"="Genus.x", "dada2_Genus"="Genus.y",
                "paratha_Species"="Species.x", "dada2_Species"="Species.y")

topBugsDada2 <- dat1 %>% group_by(dada2_Genus) %>% dplyr::summarize(count=sum(SRR3225703.y)) 
topBugsDada2 <- topBugsDada2 %>% arrange(desc(count)) 
topBugsDada2Names <- data.frame(topBugsDada2[1:30,])$dada2_Genus


dat2 <- dat1 %>% 
  group_by(paratha_Genus, dada2_Genus) %>%
  replace_na(list(paratha_Genus = "Unknown", dada2_Genus = "Unknown")) %>%
  dplyr::count() %>% 
  filter(dada2_Genus %in% topBugsDada2Names) 

dat3 <- dat2 %>%
  pivot_wider(names_from = paratha_Genus, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "dada2_Genus")

pheatmap(dat3, scale="row")

topBugsDada2 <- dat1 %>% group_by(dada2_Species) %>% dplyr::summarize(count=sum(SRR3225703.y)) 
topBugsDada2 <- topBugsDada2 %>% arrange(desc(count)) 
topBugsDada2Names <- data.frame(topBugsDada2[1:30,])$dada2_Species


dat2 <- dat1 %>% 
  group_by(paratha_Species, dada2_Species) %>%
  replace_na(list(paratha_Species = "Unknown", dada2_Species = "Unknown")) %>%
  dplyr::count() %>% 
  filter(dada2_Species %in% topBugsDada2Names) 

dat3 <- dat2 %>%
  pivot_wider(names_from = paratha_Species, values_from = n) %>%
  replace(is.na(.), 0) %>%
  column_to_rownames(var = "dada2_Species")
pheatmap(dat3, scale="row")




####################################
# Create dada2 compatible database
# Start with degapped seed db
####################################

seed.db <- seqinr::read.fasta("~/proj/parathaa/input/silva.seed_v138_1.ng.fasta", forceDNAtolower = FALSE,
                      as.string=TRUE)
#Extract names: "Accession_number.Arb_ID"
nam <- getName(seed.db)
# Extract Accession numbers
nam2 <- strsplit(nam, split="[.]")
nam3 <- sapply(nam2,"[[",1)

## Get species assignment file for larger SILVA db
silva.sp <- seqinr::read.fasta("/Users/mis696/proj/PPITAA/input/silva_species_assignment_v138.1.fa.gz" ,forceDNAtolower = FALSE,
                       as.string=TRUE)
## Extract names: "Accession.number.start.end"
silva.nam <- getName(silva.sp)
silva.nam2 <- strsplit(silva.nam, split="[.]")
# Extract Accession numbers
silva.nam3 <- sapply(silva.nam2,"[[",1)

# How many seed database sequences are found in the silva species naming file?
length(which(nam3 %in% silva.nam3))

# Which seqs are in seed db but not in silva_species_assignment training file?
excluded <- seed.db[which(!nam3 %in% silva.nam3)]

## Short answer: Eukaryotes, and species without names (some of which have numbers, i.e. Paenibacillus sp. Cp_S316	45277)
## We want some of these (the named ones that are not Eukaryotes), so we try a different approach

# Read in taxonomy file from full SILVA db
##tax <- read.delim("input/silva.seed_v138_1.tax", header=FALSE, row.names = 1, col.names = c("Name", "Taxonomy"))

inFileTaxdata <- "/Users/mis696/proj/16s-region-checker/input/taxmap_slv_ssu_ref_138.1.txt"


taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
taxdata <- taxdata %>%
  separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  dplyr::rename(Species = organism_name) %>%
  filter(Kingdom=="Bacteria") 

taxdata.sp <- taxdata %>% group_by(primaryAccession, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% tally()

#Subset to primaryAccessions from seed db
taxmatch <- taxdata.sp[which(taxdata.sp$primaryAccession %in% nam3),] 
taxmatch.df <- as.data.frame(taxmatch)
rownames(taxmatch.df) <- taxmatch.df$primaryAccession
taxmatch2 <- taxmatch.df[nam3,]

#Out of curiosity, find doubles in taxdata
#doubles <- names(table(taxdata.sp$primaryAccession)[which(table(taxdata.sp$primaryAccession)>1)])
#taxdata.sp[which(taxdata.sp$primaryAccession %in% doubles),] %>% View()

## More seed sequences show up in this full database than in the pre-trained database (based on NR99)
dim(taxmatch2)
dim(taxmatch2[!is.na(taxmatch2$Species),])

## Write seed db in training file format for DADA2 --- FIX BELOW, LENGTHS DONT MATCH
seed.db.species <- seed.db[which(!is.na(taxmatch2$Species))]
seed.db.species.names <- paste(nam3[which(!is.na(taxmatch2$Species))], taxmatch2$Species[which(!is.na(taxmatch2$Species))])

write.fasta(sequences = seed.db.species, names=seed.db.species.names, nbchar=80, file.out = "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")

### Get full taxonomies for seed database (Bacteria only)
bacRows <- which(!is.na(taxmatch2$Kingdom))
seed.db.bac <- seed.db[bacRows]
seed.db.bac.names <- apply(cbind(taxmatch2$Kingdom[bacRows],
                   taxmatch2$Phylum[bacRows],
                   taxmatch2$Class[bacRows],
                   taxmatch2$Order[bacRows],
                   taxmatch2$Family[bacRows],
                   taxmatch2$Genus[bacRows]), 1, 
             function(x) paste(x[!is.na(x)], collapse = ";"))


write.fasta(sequences = seed.db.bac , names = seed.db.bac.names, nbchar = 80, file.out = "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta")


taxa <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SRR3225703.fasta", 
                       "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta",
                       multithread=TRUE)
nChars <- grep("N", rownames(taxa))
print(paste("Removing", length(nChars), "sequences with N bases"))
taxa <- taxa[-nChars,]
taxa.sp <- addSpecies(taxa, "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")


