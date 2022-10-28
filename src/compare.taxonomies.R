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

## Read in data with taxonomy from Kelsey's analysis with SILVA
kelseyDataRaw <- read.delim("~/proj/PPITAA/input/all_samples_taxonomy_closed_reference.tsv", 
                             sep='\t', fill=T, stringsAsFactors = F, header=T)
tax_oldDADA2 <- kelseyDataRaw %>%
  select(taxonomy) %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family",
                            "Genus", "Species"), sep=";")

## Remove samples as instructed by Kelsey
kelseyDataFull <- kelseyDataRaw %>%
  select(-c(X,
            CBSF036_S66_L001, 
            CBSF060_S82_L001, 
            CBSF078_S77_L001, 
            CBSF079_S65_L001,
            Undetermined_S0_L001,
            taxonomy))

taxa <- assignTaxonomy("~/proj/PPITAA/input/amplicons_fromKelsey.fasta", 
                       "/Users/mis696/proj/PPITAA/input/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=TRUE)
taxa.sp <- addSpecies(taxa, "/Users/mis696/proj/PPITAA/input/silva_species_assignment_v138.1.fa.gz")


tax_dada <- as.matrix(taxa.sp)
#Get IDs from format_amplicon_table.R "kelseyData" object
rownames(tax_dada) <- paste0(kelseyData$taxaIDs, "a")

TAX_dada <- tax_table(tax_dada)

otumat <- as.matrix(kelseyDataFull)
rownames(otumat) <- paste0(kelseyData$taxaIDs, "a")
colnames(otumat) <- paste0(colnames(otumat), "a")
OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)

samp_dada <- data.frame(colnames(OTU), rep("DADA2", length(colnames(OTU))))

colnames(samp_dada) <- c("sampleID", "Taxonomy_type")
rownames(samp_dada) <- paste0(samp_dada$sampleID, "a")
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
ps1.com <- ps1_parathaa

# We need to set Palette
taxic <- as.data.frame(ps1.com@tax_table) # this will help in setting large color options

#getPalette = colorRampPalette(brewer.pal(colourCount, ))  # change the palette as well as the number of colors will change according to palette.

taxic$OTU <- rownames(taxic) # Add the OTU ids from OTU table into the taxa table at the end.
colnames(taxic) # You can see that we now have extra taxonomy levels.

taxmat <- as.matrix(taxic) # convert it into a matrix.
new.tax <- tax_table(taxmat) # convert into phyloseq compatible file.
tax_table(ps1.com) <- new.tax # incroporate into phyloseq Object


# now edit the unclassified taxa
tax_table(ps1.com)[tax_table(ps1.com)[, "Family"] == "", "Family"] <- "Unclassified family"
tax_table(ps1.com)[tax_table(ps1.com)[, "Phylum"] == "", "Phylum"] <- "Unclassified phylum"

# it would be nice to have the Taxonomic names in italics.
# for that we set this

guide_italics <- guides(fill = guide_legend(label.theme = element_text(
  size = 15,
  face = "italic", colour = "Black", angle = 0
)))


## Now we need to plot at family level, we can do it as follows:

# first remove the phy_tree

ps1.com@phy_tree <- NULL


Levels <- c("Phylum", "Class", "Order", "Family", "Genus")

# Taxonomy plots
for(level in Levels){
  ps1.com.lev <- aggregate_rare(ps1.com, level, detection = .01/100, prevalence = 10/100)
  ps1.com.lev.agg <- aggregate_top_taxa2(ps1.com.lev, 10, level)
  ps1.com.lev.agg.rel <- microbiome::transform(ps1.com.lev.agg, "compositional")
  plot.composition.relAbun <- plot_composition(ps1.com.lev.agg.rel,
                                               sample.sort = "sampleID",
                                               #sample.sort = "Taxonomy_type",
                                               otu.sort = "abundance",
                                               x.label = "sampleID") 
  plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
  plot.composition.relAbun <- plot.composition.relAbun +  theme_bw()  +scale_fill_brewer(level, palette = "Paired") #+ theme_bw() 
  plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
  plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + guide_italics + theme(legend.title = element_text(size = 18))
  ggsave(filename=paste0("~/proj/PPITAA/output/viz/", level, "parathaa_taxonomy.png"), plot.composition.relAbun, width=10, height=6)
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

ggsave(filename=paste0("~/proj/PPITAA/output/viz/AllLevels_PCoA_bray.png"), plot=p1, height = 8, width=6, dpi=300)



#Heatmap
dat1 <- bind_cols(data.frame(tax_table(ps1_parathaa))$Genus, data.frame(tax_table(ps1_dada))$Genus)

topBugsDada2 <- dat1 %>% group_by(dada2_Genus) %>% count() 
topBugsDada2 <- topBugsDada2 %>% arrange(desc(n)) 
topBugsDada2Names <- data.frame(topBugsDada2[1:30,])$dada2_Genus


colnames(dat1) <- c("paratha_Genus", "dada2_Genus")
dat2 <- dat1 %>% 
  group_by(paratha_Genus, dada2_Genus) %>%
  replace_na(list(paratha_Genus = "Unknown", dada2_Genus = "Unknown")) %>%
  count() %>% 
  filter(dada2_Genus %in% topBugsDada2Names) 

  

  
