## Compare taxonomies from SILVA vs PARATHAA

library(microbiome) # data analysis and visualisation
library(phyloseq) # also the basis of data object. Data analysis and visualisation
library(RColorBrewer) # nice color options
library(ggpubr) # publication quality figures, based on ggplot2
library(tidyverse) # data handling  
library(treeio)
library(gridExtra)
library(lemon)
library(vegan)
library(microbiomeutilities)
library(pheatmap)
library(dada2)
library(devtools)
library(pals)
library(ComplexHeatmap)
suppressPackageStartupMessages(library(seqinr))
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
#taxaNR99 <- assignTaxonomy("~/proj/PPITAA/input/amplicons_fromKelsey.fasta", 
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

samp_dada <- data.frame(colnames(otumat), rep("DADA2", length(colnames(otumat))))

colnames(otumat) <- paste0(colnames(otumat), "a")
OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)


colnames(samp_dada) <- c("sampleID", "Taxonomy_type")
rownames(samp_dada) <-  paste0(samp_dada$sampleID, "a")
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


ps1_all<- merge_phyloseq(ps1_dada, ps1_parathaa, ps1_spingo)
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
tax_table(ps1.com)[is.na(tax_table(ps1.com)[, "Species"]), "Species"] <- "Unknown"
tax_table(ps1.com)[which(tax_table(ps1.com)[, "Species"]=="UNCLASSIFIED"), "Species"] <- "Unknown"
tax_table(ps1.com)[which(tax_table(ps1.com)[, "Species"]=="AMBIGUOUS"), "Species"] <- "Ambiguous"


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
Levels <- "Species"

# Taxonomy plots
for(level in Levels){

  ps1.com.rel <- microbiome::transform(ps1.com, "compositional")
  ps1.com.rel.lev <- aggregate_rare(ps1.com.rel, level, detection = .1/100, prevalence = 10/100)
  ps1.com.rel.lev.agg <- aggregate_top_taxa2(ps1.com.rel.lev, 15, level)

  
  #ps1.com.rel.lev.agg <- subset_samples(ps1.com.rel.lev.agg, Taxonomy_type!="SPINGO")
  
  #ps1.com.lev <- ps1.com

  plot.composition.relAbun <- plot_composition(ps1.com.rel.lev.agg,
                                               #sample.sort = "sampleID",
                                               sample.sort = "Region",
                                               #otu.sort = "abundance",
                                               x.label = "sampleID") 
  #xlabs <- rep("", nrow(sample_data(ps1.com.rel.lev.agg)))
  #xlabs[ceiling(nrow(sample_data(ps1.com.rel.lev.agg))/4)] <- "DADA2"
  #xlabs[ceiling(nrow(sample_data(ps1.com.rel.lev.agg))*3/4)+1] <- "parathaa"
  #xlabs <- c("DADA2", "Parathaa")
  xlabs <- c("DADA2 V1V2", "Parathaa V1V2", "SPINGO V1V2", "DADA2 V4V5", "Parathaa V4V5", "SPINGO V4V5")
  plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom") 
  plot.composition.relAbun <- plot.composition.relAbun +  theme_bw()  + scale_fill_manual(level, values = stepped2())# scale_fill_brewer(palette="Paired")  #
  #plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90)) 
  plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + #guide_italics + theme(legend.title = element_text(size = 18)) + 
    guides(fill= 
             guide_legend(
               ncol=1,
               label.theme = element_text(face="italic", size=14),
               title.theme = element_text(size=18),
               )
           )
  plot.composition.relAbun <- plot.composition.relAbun +   theme(axis.text.x = element_text(size=12, angle=45, hjust=1)) + scale_x_discrete(labels=xlabs) 
  #ggsave(filename=paste0("~/proj/parathaa/output/viz/", level, "parathaa_taxonomy.png"), plot.composition.relAbun, width=10, height=6)
  if(level!="Species")
    ggsave(filename=paste0("/Users/mis696/proj/parathaa/output/20221231_Mock", level, "taxonomy.png"), plot.composition.relAbun, width=6, height=6)
  if(level=="Species")
    ggsave(filename=paste0("/Users/mis696/proj/parathaa/output/20221231_Mock", level, "taxonomy.png"), plot.composition.relAbun, width=10, height=6)
}

## Heatmap of mock community
ps1.com.rel <- microbiome::transform(ps1.com, "compositional")
ps1.com.rel.lev <- aggregate_rare(ps1.com.rel, level, detection = .1/100, prevalence = 10/100)
ps1.com.rel.lev.agg <- aggregate_taxa(ps1.com.rel.lev, level)
forHeat <- t(as.matrix(otu_table(ps1.com.rel.lev.agg)))
rownames(forHeat) <- c("V1V2 DADA2", "V4V5 DADA2", 
                       "V1V2 Parathaa", "V4V5 Parathaa", 
                       "V1V2 SPINGO", "V4V5 SPINGO")
forHeat <- forHeat[c(1,3,5,2,4,6),]

#View(forHeat)
included <- c("Actinomyces odontolyticus",
              "Bacillus cereus",
              "Bacillus cereus;Bacillus thuringiensis",
              "Bacillus anthracis;Bacillus cereus;Bacillus thuringiensis",
              "Bacteroides vulgatus"                                     ,
              "Clostridium diolis"                                        ,
              "Clostridium sensu stricto 1 diolis"                        ,
              "Clostridium beijerinckii",
              "Cutibacterium acnes",
              "Deinococcus radiodurans"                                  ,
              "Enterococcus canis;Enterococcus faecalis"         ,
              "Enterococcus faecalis"                                    ,
              "Lactobacillus gasseri;Lactobacillus johnsonii"            ,
              "Lactobacillus gasseri",
              "Listeria innocua;Listeria monocytogenes"                   ,
              "Listeria monocytogenes"                                   ,
              "Neisseria meningitidis"                                   ,
              "Propionibacterium acnes", 
              "Pseudomonas aeruginosa"                                  ,
              "Staphylococcus aureus"                                    ,
              "Staphylococcus epidermidis"  ,
              "Staphylococcus aureus;Staphylococcus epidermidis"         ,  
              "Streptococcus agalactiae"                                 ,
              "Streptococcus mutans"                                     ,
              "Streptococcus pneumoniae"    )
others <- c("Ambiguous", "Unknown", "Other")
included.df <- data.frame(rep("Not included", length(colnames(forHeat))))
colnames(included.df) <- "Mock Community"
included.df[which(colnames(forHeat) %in% included), "Mock Community"] <- "Included"
included.df[which(colnames(forHeat) %in% others), "Mock Community"] <- "Other"
rownames(included.df) <- colnames(forHeat)

type <- gsub("s\\d+_", "", colnames(forHeat))
col_fun = viridis(3)
ha = columnAnnotation(
  df = included.df[,1], show_annotation_name=FALSE, 
  annotation_legend_param = list(title = "Mock Community"),
  col=list(df = c("Included" =  viridis(4)[2],"Not included" = viridis(4)[3], "Other" = viridis(4)[1]))
)
region.df <- data.frame(c(rep("V1V2", length(rownames(forHeat))/2), rep("V4V5", length(rownames(forHeat))/2)))
colnames(region.df) <- "16S Region"
rownames(region.df) <- rownames(forHeat)
  
ra = rowAnnotation(
  df = region.df[,1], show_annotation_name=FALSE, 
  annotation_legend_param = list(title = "16S Region"),
  col=list(df = c("V1V2" =  turbo(7)[2],"V4V5" = turbo(7)[4]))
)
png("~/proj/parathaa/output/viz/Fig_2A_MockHeatmap_20221228.png",width=10,height=5.5,units="in",res=600)
ht1 <- Heatmap(sqrt(forHeat), name = "sqrt(Rel.\nabundance)", col=c("grey", rev(magma(100))),
               column_names_max_height = unit(10, "cm"),
               cluster_columns = FALSE, cluster_rows=FALSE, column_names_rot = 45, 
               column_names_side = "top", column_names_gp = grid::gpar(fontsize = 10), 
               column_split = included.df[,"Mock Community"], column_title=NULL, top_annotation = ha,
               row_split = region.df[,"16S Region"], row_title=NULL, right_annotation = ra,
               row_names_gp = grid::gpar(fontsize = 10))  %v% NULL

draw(ht1)

dev.off()

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
    geom_polygon(aes(group=sampleID), color="grey") + theme_bw()
    
  plotList[[level]] <- p1
}

p1 <- ggarrange(plotList[["Phylum"]]+theme(legend.position='hidden'), plotList[["Class"]]+theme(legend.position='hidden'), 
          plotList[["Order"]]+theme(legend.position='hidden'), plotList[["Family"]]+theme(legend.position='hidden'),
          plotList[["Genus"]]+theme(legend.position='hidden'), plotList[["Species"]], ncol=2, nrow=3, common.legend = TRUE, legend="bottom") 

ggsave(filename=paste0("~/proj/parathaa/output/viz/20221214_Mockv4v5_PCoA_bray_seedDada.png"), plot=p1, height = 6, width=3, dpi=300)


#Heatmap
dat1 <- bind_cols(data.frame(tax_table(ps1_parathaa.V4V5))$Genus, data.frame(tax_table(ps1_dada.V4V5))$Genus)
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

pheatmap(as.matrix(dat3))#, scale="row")


###############################
### Mock community analysis ###
###############################

## Assign taxonomy with dada2 and full SILVA db
if(FALSE){
taxaV4V5 <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SRR3225703.fasta", 
                       "/Users/mis696/proj/PPITAA/input/silva_nr99_v138.1_train_set.fa.gz",
                       multithread=TRUE)
nChars <- grep("N", rownames(taxaV4V5))
print(paste("Removing", length(nChars), "sequences with N bases"))
taxaV4V5 <- taxaV4V5[-nChars,]
taxa.sp <- addSpecies(taxa, "/Users/mis696/proj/PPITAA/input/silva_species_assignment_v138.1.fa.gz")
}

## Assign taxonomy with dada2 and seed SILVA db
## V4V5
taxaV4V5 <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SRR3225703.fasta", 
                           "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta",
                       multithread=TRUE)
nChars <- grep("N", rownames(taxaV4V5))
print(paste("Removing", length(nChars), "sequences with N bases"))
taxaV4V5 <- taxaV4V5[-nChars,]
taxaV4V5.sp <- addSpecies(taxaV4V5, "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")

otutabV4V5 <- cbind(rownames(taxaV4V5.sp), taxaV4V5.sp)
otutabV4V5 <- as_tibble(otutabV4V5) %>% dplyr::rename("ASV"="V1") %>% dplyr::count(ASV, name="SRR3225703.V4V5") 
otumatV4V5 <- as.matrix(otutabV4V5$SRR3225703.V4V5)
rownames(otumatV4V5) <- otutabV4V5$ASV
colnames(otumatV4V5) <- "SRR3225703.V4V5"

OTU_dada.V4V5 <- otu_table(otumatV4V5, taxa_are_rows = TRUE)

taxtabV4V5 <- cbind(rownames(taxaV4V5.sp), taxaV4V5.sp)
taxtabV4V5 <- as_tibble(taxtabV4V5) %>% 
  dplyr::rename("ASV"="V1") %>% 
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225703.V4V5") 
taxmatV4V5 <- taxtabV4V5 %>% as.matrix
rownames(taxmatV4V5) <- taxtabV4V5$ASV
taxmatV4V5 <- taxmatV4V5[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
taxmatV4V5 <- cbind(taxmatV4V5, "Species2"=NA)
taxmatV4V5[which(!is.na(taxmatV4V5[,"Species"])), "Species2"] <-  
  paste(taxmatV4V5[which(!is.na(taxmatV4V5[,"Species"])), "Genus"], 
        taxmatV4V5[which(!is.na(taxmatV4V5[,"Species"])), "Species"] )
taxmatV4V5 <- taxmatV4V5[,!colnames(taxmatV4V5) %in% "Species"]
colnames(taxmatV4V5)[which(colnames(taxmatV4V5)=="Species2")] <- "Species"

TAX_dada.V4V5 <- tax_table(taxmatV4V5)

samp_dada.V4V5 <- data.frame(colnames(OTU_dada.V4V5), 
                             rep("DADA2", length(colnames(OTU_dada.V4V5))),
                             rep("V4V5", length(colnames(OTU_dada.V4V5))))
colnames(samp_dada.V4V5) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_dada.V4V5) <- samp_dada.V4V5$sampleID
str(samp_dada.V4V5)
SAMP_dada.V4V5 <- sample_data(samp_dada.V4V5)

ps1_dada.V4V5 <- phyloseq(OTU_dada.V4V5, TAX_dada.V4V5, SAMP_dada.V4V5)

## V1V2 -- in two parts for memory reasons
taxaV1V2.0 <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SRR3225701.0.fasta", 
                           "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta",
                           multithread=TRUE)
nChars <- grep("N", rownames(taxaV1V2.0))
print(paste("Removing", length(nChars), "sequences with N bases"))
#taxaV1V2.0 <- taxaV1V2.0[-nChars,]
taxaV1V2.0.sp <- addSpecies(taxaV1V2.0, "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")

taxaV1V2.1 <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SRR3225701.1.fasta", 
                             "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta",
                             multithread=TRUE)
taxaV1V2.1.sp <- addSpecies(taxaV1V2.1, "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")

otutabV1V2.0 <- cbind(rownames(taxaV1V2.0.sp), taxaV1V2.0.sp)
otutabV1V2.1 <- cbind(rownames(taxaV1V2.1.sp), taxaV1V2.1.sp)
otutabV1V2 <- bind_rows(as_tibble(otutabV1V2.0), as_tibble(otutabV1V2.1))
otutabV1V2 <- otutabV1V2 %>% dplyr::rename("ASV"="V1") %>% dplyr::count(ASV, name="SRR3225701.V1V2") 
otumatV1V2 <- as.matrix(otutabV1V2$SRR3225701.V1V2)
rownames(otumatV1V2) <- otutabV1V2$ASV
colnames(otumatV1V2) <- "SRR3225701.V1V2"

OTU_dada.V1V2 <- otu_table(otumatV1V2, taxa_are_rows = TRUE)

taxtabV1V2 <- rbind(
  cbind(rownames(taxaV1V2.0.sp), taxaV1V2.0.sp),
  cbind(rownames(taxaV1V2.1.sp), taxaV1V2.1.sp)
)
taxtabV1V2 <- as_tibble(taxtabV1V2) %>% 
  dplyr::rename("ASV"="V1") %>% 
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225701.V1V2")
taxtabV1V2 <- taxtabV1V2 %>%
  group_by(ASV) %>% slice_max(SRR3225701.V1V2, with_ties = FALSE)

taxmatV1V2 <- taxtabV1V2 %>% as.matrix
rownames(taxmatV1V2) <- taxtabV1V2$ASV
taxmatV1V2 <- taxmatV1V2[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
taxmatV1V2 <- cbind(taxmatV1V2, "Species2"=NA)
taxmatV1V2[which(!is.na(taxmatV1V2[,"Species"])), "Species2"] <-  
  paste(taxmatV1V2[which(!is.na(taxmatV1V2[,"Species"])), "Genus"], 
        taxmatV1V2[which(!is.na(taxmatV1V2[,"Species"])), "Species"] )
taxmatV1V2 <- taxmatV1V2[,!colnames(taxmatV1V2) %in% "Species"]
colnames(taxmatV1V2)[which(colnames(taxmatV1V2)=="Species2")] <- "Species"

TAX_dada.V1V2 <- tax_table(taxmatV1V2)

samp_dada.V1V2 <- data.frame(colnames(OTU_dada.V1V2), 
                             rep("DADA2", length(colnames(OTU_dada.V1V2))),
                             rep("V1V2", length(colnames(OTU_dada.V1V2))))
colnames(samp_dada.V1V2) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_dada.V1V2) <- samp_dada.V1V2$sampleID
str(samp_dada.V1V2)
SAMP_dada.V1V2 <- sample_data(samp_dada.V1V2)

ps1_dada.V1V2 <- phyloseq(OTU_dada.V1V2, TAX_dada.V1V2, SAMP_dada.V1V2)




## Read in parathaa data
#parathaData <- read.delim("/Users/mis696/proj/parathaa/output/20221101_MiSeqV4V5Mock/taxonomic_assignments.tsv", 
#                          sep='\t', fill=T, stringsAsFactors = F, header=T)
parathaData.V4V5 <- read.delim("/Users/mis696/proj/parathaa/output/20221215_MiSeqV4V5Mock_Sp/taxonomic_assignments.tsv", 
                          sep='\t', fill=T, stringsAsFactors = F, header=T)

tax_paratha.V4V5 <- parathaData.V4V5 %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  group_by(query.name) %>% 
  filter(row_number()==1)

## utility to refine species names:
refine_species_names <-function(x){
  x <- x %>% mutate(Species = recode(Species,
                                     "Bacillus anthracis;Bacillus cereus;Bacillus cereus ATCC 10987;Bacillus thuringiensis" = "Bacillus anthracis;Bacillus cereus;Bacillus thuringiensis",
                                     "Enterococcus canis;Enterococcus faecalis V583;uncultured bacterium" = "Enterococcus canis;Enterococcus faecalis;unknown",
                                     "Lactobacillus gasseri ATCC 33323 = JCM 1131;Lactobacillus johnsonii N6.2" = "Lactobacillus gasseri;Lactobacillus johnsonii",
                                     "Lactobacillus gasseri;Lactobacillus johnsonii N6.2" = "Lactobacillus gasseri;Lactobacillus johnsonii",
                                     "Listeria innocua;Listeria monocytogenes;Listeria monocytogenes N53-1" = "Listeria innocua;Listeria monocytogenes",
                                     "Neisseria meningitidis alpha14;Neisseria meningitidis MC58;Neisseria meningitidis Z2491" = "Neisseria meningitidis",
                                     "Streptococcus agalactiae 18RS21" = "Streptococcus agalactiae",
                                     "Streptococcus mutans;Streptococcus mutans UA159" = "Streptococcus mutans",
                                     "Streptococcus pneumoniae;Streptococcus pneumoniae TIGR4" = "Streptococcus pneumoniae"
  ))
  return(x)
}


parsed = seqinr::read.fasta(file('/Users/mis696/proj/parathaa/input/SRR3225703.fasta'), as.string = TRUE,
                    forceDNAtolower = FALSE, whole.header = FALSE)
table = data.frame("ASV"=unlist(parsed), "query.name" = sapply(parsed, attr, 'name'), row.names=NULL)
#View(tax_paratha)

tax_paratha.V4V5 <- merge(tax_paratha.V4V5, table)  %>%
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225703.V4V5") 

#tax_paratha.V4V5 <- refine_species_names(tax_paratha.V4V5)

taxmatV4V5 <- tax_paratha.V4V5 %>% as.matrix
rownames(taxmatV4V5) <- paste0("x", tax_paratha.V4V5$ASV )
taxmatV4V5 <- taxmatV4V5[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa.V4V5 <- tax_table(taxmatV4V5)

otumatV4V5 <- as.matrix(tax_paratha.V4V5$SRR3225703.V4V5)
rownames(otumatV4V5) <- paste0("x", tax_paratha.V4V5$ASV )
colnames(otumatV4V5) <- "SRR3225703b.V4V5"
OTU_parathaa.V4V5 <- otu_table(otumatV4V5, taxa_are_rows = TRUE)

samp_parathaa.V4V5 <- data.frame(colnames(OTU_parathaa.V4V5), 
                             rep("Parathaa", length(colnames(OTU_parathaa.V4V5))),
                             rep("V4V5", length(colnames(OTU_parathaa.V4V5))))
colnames(samp_parathaa.V4V5) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_parathaa.V4V5) <- samp_parathaa.V4V5$sampleID
str(samp_parathaa.V4V5)
SAMP_parathaa.V4V5 <- sample_data(samp_parathaa.V4V5)

ps1_parathaa.V4V5 <- phyloseq(OTU_parathaa.V4V5, TAX_parathaa.V4V5, SAMP_parathaa.V4V5)

## V1V2
parathaData.V1V2 <- read.delim("/Users/mis696/proj/parathaa/output/20221219_MiSeqV1V2Mock/taxonomic_assignments.tsv", 
                               sep='\t', fill=T, stringsAsFactors = F, header=T)

tax_paratha.V1V2 <- parathaData.V1V2 %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  group_by(query.name) %>% 
  filter(row_number()==1)

## utility to refine species names:
refine_species_names <-function(x){
  x <- x %>% mutate(Species = recode(Species,
                                     "Bacillus anthracis;Bacillus cereus;Bacillus cereus ATCC 10987;Bacillus thuringiensis" = "Bacillus anthracis;Bacillus cereus;Bacillus thuringiensis",
                                     "Enterococcus canis;Enterococcus faecalis V583;uncultured bacterium" = "Enterococcus canis;Enterococcus faecalis;unknown",
                                     "Lactobacillus gasseri ATCC 33323 = JCM 1131;Lactobacillus johnsonii N6.2" = "Lactobacillus gasseri;Lactobacillus johnsonii",
                                     "Lactobacillus gasseri;Lactobacillus johnsonii N6.2" = "Lactobacillus gasseri;Lactobacillus johnsonii",
                                     "Listeria innocua;Listeria monocytogenes;Listeria monocytogenes N53-1" = "Listeria innocua;Listeria monocytogenes",
                                     "Neisseria meningitidis alpha14;Neisseria meningitidis MC58;Neisseria meningitidis Z2491" = "Neisseria meningitidis",
                                     "Streptococcus agalactiae 18RS21" = "Streptococcus agalactiae",
                                     "Streptococcus mutans;Streptococcus mutans UA159" = "Streptococcus mutans",
                                     "Streptococcus pneumoniae;Streptococcus pneumoniae TIGR4" = "Streptococcus pneumoniae",
                                     "Synechococcus sp.;Synechococcus sp. PCC 7002;Synechococcus sp. PH40" = "Synechococcus sp. PCC 7002;Synechococcus sp. PH40"
  ))
  return(x)
}


parsed = seqinr::read.fasta(file('/Users/mis696/proj/parathaa/input/SRR3225701.fasta'), as.string = TRUE,
                            forceDNAtolower = FALSE, whole.header = FALSE)
table = data.frame("ASV"=unlist(parsed), "query.name" = sapply(parsed, attr, 'name'), row.names=NULL)
#View(tax_paratha)

tax_paratha.V1V2 <- merge(tax_paratha.V1V2, table)  %>%
  dplyr::group_by(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(ASV, name="SRR3225701.V1V2") 

tax_paratha.V1V2 <- refine_species_names(tax_paratha.V1V2)

taxmatV1V2 <- tax_paratha.V1V2 %>% as.matrix
rownames(taxmatV1V2) <- paste0("x", tax_paratha.V1V2$ASV )
taxmatV1V2 <- taxmatV1V2[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa.V1V2 <- tax_table(taxmatV1V2)

otumatV1V2 <- as.matrix(tax_paratha.V1V2$SRR3225701.V1V2)
rownames(otumatV1V2) <- paste0("x", tax_paratha.V1V2$ASV )
colnames(otumatV1V2) <- "SRR3225701b.V1V2"
OTU_parathaa.V1V2 <- otu_table(otumatV1V2, taxa_are_rows = TRUE)

samp_parathaa.V1V2 <- data.frame(colnames(OTU_parathaa.V1V2), 
                                 rep("Parathaa", length(colnames(OTU_parathaa.V1V2))),
                                 rep("V1V2", length(colnames(OTU_parathaa.V1V2))))
colnames(samp_parathaa.V1V2) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_parathaa.V1V2) <- samp_parathaa.V1V2$sampleID
str(samp_parathaa.V1V2)
SAMP_parathaa.V1V2 <- sample_data(samp_parathaa.V1V2)

ps1_parathaa.V1V2 <- phyloseq(OTU_parathaa.V1V2, TAX_parathaa.V1V2, SAMP_parathaa.V1V2)



## Read in SPINGO results
spingo.V4V5 <- read.delim("/Users/mis696/proj/parathaa/output/20221101_MiSeqV4V5Mock/SRR3225703.summary.txt",
                     header = F, col.names = c("Species", "SRR3225703c"))
spingo.V4V5$Species <- str_replace_all(spingo.V4V5$Species, "\\(", "")
spingo.V4V5$SRR3225703c <- as.numeric(str_replace_all(spingo.V4V5$SRR3225703c, "\\)", ""))

taxmat <- as.matrix(spingo.V4V5$Species)
rownames(taxmat) <- spingo.V4V5$Species 
colnames(taxmat) <- "Species"
taxmat[,"Species"] <- str_replace_all(taxmat[,"Species"], "_", " ")
TAX_spingo <- tax_table(taxmat)

otumat <- as.matrix(spingo.V4V5$SRR3225703c)
rownames(otumat) <- spingo.V4V5$Species
colnames(otumat) <- "SRR3225703c"
OTU_spingo <- otu_table(otumat, taxa_are_rows = TRUE)

samp_spingo <- data.frame(colnames(OTU_spingo), rep("SPINGO", length(colnames(OTU_spingo))),
                          rep("V4V5", length(colnames(OTU_spingo))))
colnames(samp_spingo) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_spingo) <- samp_spingo$sampleID
str(samp_spingo)
SAMP_spingo <- sample_data(samp_spingo)

ps1_spingo.V4V5 <- phyloseq(OTU_spingo, TAX_spingo, SAMP_spingo)

spingo.V1V2 <- read.delim("/Users/mis696/proj/parathaa/output/20221219_MiSeqV1V2Mock/SRR3225701.summary.txt",
                          header = F, col.names = c("Species", "SRR3225701c"))
spingo.V1V2$Species <- str_replace_all(spingo.V1V2$Species, "\\(", "")
spingo.V1V2$SRR3225701c <- as.numeric(str_replace_all(spingo.V1V2$SRR3225701c, "\\)", ""))

taxmat <- as.matrix(spingo.V1V2$Species)
rownames(taxmat) <- spingo.V1V2$Species 
colnames(taxmat) <- "Species"
taxmat[,"Species"] <- str_replace_all(taxmat[,"Species"], "_", " ")
TAX_spingo <- tax_table(taxmat)

otumat <- as.matrix(spingo.V1V2$SRR3225701c)
rownames(otumat) <- spingo.V1V2$Species
colnames(otumat) <- "SRR3225701c"
OTU_spingo <- otu_table(otumat, taxa_are_rows = TRUE)

samp_spingo <- data.frame(colnames(OTU_spingo), rep("SPINGO", length(colnames(OTU_spingo))),
                          rep("V1V2", length(colnames(OTU_spingo))))
colnames(samp_spingo) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_spingo) <- samp_spingo$sampleID
str(samp_spingo)
SAMP_spingo <- sample_data(samp_spingo)

ps1_spingo.V1V2 <- phyloseq(OTU_spingo, TAX_spingo, SAMP_spingo)

ps1_all<- merge_phyloseq(ps1_dada.V1V2, ps1_dada.V4V5, 
                         ps1_parathaa.V1V2, ps1_parathaa.V4V5, 
                         ps1_spingo.V1V2, ps1_spingo.V4V5)


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

## Read in seed db
seed.db <- seqinr::read.fasta("~/proj/parathaa/input/silva.seed_v138_1.ng.fasta", forceDNAtolower = FALSE,
                      as.string=TRUE)
#Extract names: "Accession_number.Arb_ID"
nam <- getName(seed.db)
# Extract Accession numbers
nam2 <- strsplit(nam, split="[.]")
nam3 <- sapply(nam2,"[[",1)

## Get dada2 species assignment file for larger SILVA db
silva.sp <- seqinr::read.fasta("/Users/mis696/proj/PPITAA/input/silva_species_assignment_v138.1.fa.gz" ,forceDNAtolower = FALSE,
                       as.string=TRUE)
## Extract names: "Accession.number.start.end"
silva.nam <- getName(silva.sp)
silva.nam2 <- strsplit(silva.nam, split="[.]")
# Extract Accession numbers
silva.nam3 <- sapply(silva.nam2,"[[",1)

# How many seed database sequences are found in the dada2 silva species naming file?
length(which(nam3 %in% silva.nam3))

# Which seqs are in seed db but not in silva_species_assignment training file?
excluded <- seed.db[which(!nam3 %in% silva.nam3)]
excluded <- nam3[!nam3 %in% silva.nam3]

## Short answer: Eukaryotes, and species without names (some of which have numbers, i.e. Paenibacillus sp. Cp_S316	45277)
## We want some of these (the named ones that are not Eukaryotes), so we try a different approach

### Before doing this, going to try using only those from sp naming file:
silva.sp.seed <- silva.sp[silva.nam3 %in% nam3]
silva.sp.seed.names <- substring(getAnnot(silva.sp.seed), 2)

write.fasta(silva.sp.seed, names=silva.sp.seed.names, nbchar=80, file.out = "/Users/mis696/proj/parathaa/input/silva_species_assignment_v138.1.seedonly.fa")

### And similarly, take full training set and subset to seed db:
silva.train <- seqinr::read.fasta("/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.fa" ,forceDNAtolower = FALSE,
                               as.string=TRUE)

silva.train.seed <- silva.train[silva.train %in% silva.sp.seed]
write.fasta(silva.train.seed, names=getName(silva.train.seed), nbchar=80, file.out = "/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.seedonly.fa")


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
taxdata <- taxdata %>%
  mutate_at(vars("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), na_if, "uncultured")
taxdata <- taxdata %>%
  mutate_at(vars("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), na_if, "")

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

## Write seed db in training file format for DADA2 
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
seed.db.bac.names <- paste0(seed.db.bac.names, ";")


write.fasta(sequences = seed.db.bac , names = seed.db.bac.names, nbchar = 80, file.out = "/Users/mis696/proj/parathaa/input/20230111.silva.seed_v138_1.ng.dada.fasta")

##################################
## Synthetic community analysis ##
##################################

# Read in parathaa data
## Read in V1V2 parathaa data
#"/Users/mis696/proj/parathaa/output/20230110_SyntheticV1V2_nameHarmonizing/taxonomic_assignments.tsv",
parathaData <- read.delim("/Users/mis696/proj/parathaa/output/20230119_SyntheticV1V2_SpThreshold003/taxonomic_assignments.tsv",
  #"/Users/mis696/proj/parathaa/output/20230120_SyntheticV1V2_removeUndefSp/taxonomic_assignments.tsv", 
                          sep='\t', fill=T, stringsAsFactors = F, header=T)
tax_paratha <- parathaData %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  group_by(query.name) %>% 
  filter(row_number()==1)

taxmat <- tax_paratha %>% as.matrix
rownames(taxmat) <- tax_paratha$query.name
taxmat <- taxmat[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa <- tax_table(taxmat)

otutab <- tax_paratha %>% 
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

ps1_parathaa.V1V2 <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)

## Read in V4V5 parathaa data
parathaData <- read.delim("/Users/mis696/proj/parathaa/output/20230109_SyntheticV4V5_nameHarmonizing/taxonomic_assignments.tsv", 
                          sep='\t', fill=T, stringsAsFactors = F, header=T)
tax_paratha <- parathaData %>%
  dplyr::select(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  group_by(query.name) %>% 
  filter(row_number()==1)

taxmat <- tax_paratha %>% as.matrix
rownames(taxmat) <- tax_paratha$query.name
taxmat <- taxmat[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]  
TAX_parathaa <- tax_table(taxmat)

otutab <- tax_paratha %>% 
  dplyr::group_by(query.name, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(query.name, name="Parathaa") 
otumat <- as.matrix(otutab$Parathaa)
rownames(otumat) <- otutab$query.name
colnames(otumat) <- "Parathaa.V4V5"
OTU_parathaa <- otu_table(otumat, taxa_are_rows = TRUE)

samp_parathaa <- data.frame(colnames(OTU_parathaa), rep("parathaa", length(colnames(OTU_parathaa))),
                            rep("V4V5", length(colnames(OTU_parathaa))))
colnames(samp_parathaa) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_parathaa) <- samp_parathaa$sampleID
str(samp_parathaa)
SAMP_parathaa <- sample_data(samp_parathaa)

ps1_parathaa.V4V5 <- phyloseq(OTU_parathaa, TAX_parathaa, SAMP_parathaa)



## Assign taxonomy with dada2
## V4V5
## First, get names and sequences from fasta file
getNames <- read.fasta(file("/Users/mis696/proj/parathaa/input/SILVAsubsample_SeedGenera_V4V5.pcr.fasta"), as.string = TRUE,
                       forceDNAtolower = FALSE, whole.header = FALSE)
names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
names1 <- names1[,1] %>%
  str_remove(">")
name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)

## Next, assign taxonomy to genus level with DADA2
set.seed(3874)
taxa <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SILVAsubsample_SeedGenera_V4V5.pcr.fasta", 
                       "/Users/mis696/proj/parathaa/input/20230111.silva.seed_v138_1.ng.dada.fasta",
                        #"/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.seedonly.fa",
                       #"/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.fa",
                       #"/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta", 
                       multithread=TRUE)
## Remove sequences with undefined ("N") bases, store until after species assignment
taxa.test <- as.data.frame(taxa)
taxa.test$taxaIDs <- names1
nChars <- grep("N", rownames(taxa.test))
print(paste("Removing", length(nChars), "sequences with N bases"))
withNbases <- taxa.test[nChars,]
taxa <- taxa[-nChars,]

## Perform species assignment with DADA2
taxa.sp <- addSpecies(taxa, "/Users/mis696/proj/parathaa/input/silva_species_assignment_v138.1.seedonly.fa")
#                      "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")

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


TAX_dada <- tax_table(tax_dada3)


otutab <- as.data.frame(tax_dada2) %>% 
  select(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::group_by(taxaIDs, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  dplyr::count(taxaIDs, name="DADA2") 
otumat <- as.matrix(otutab$DADA2)
rownames(otumat) <- otutab$taxaIDs
colnames(otumat) <- "DADA2.V4V5"
OTU_dada <- otu_table(otumat, taxa_are_rows = TRUE)

samp_dada <- data.frame(colnames(OTU_dada), rep("DADA2", length(colnames(OTU_dada))),
                        rep("V4V5", length(colnames(OTU_dada))))
colnames(samp_dada) <- c("sampleID", "Taxonomy_type", "Region")
rownames(samp_dada) <- samp_dada$sampleID
str(samp_dada)
SAMP_dada <- sample_data(samp_dada)

ps1_dada.V4V5 <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)

print(ps1_dada.V4V5)


#V1V2
## First, get names and sequences from fasta file
getNames <- read.fasta(file("/Users/mis696/proj/parathaa/input/SILVAsubsample_SeedGenera_V1V2.pcr.fasta"), as.string = TRUE,
                       forceDNAtolower = FALSE, whole.header = FALSE)
names1 <- str_split(getName(getNames), "\t", simplify=TRUE)
names1 <- names1[,1] %>%
  str_remove(">")
name.df <- data.frame("sequence" = unlist(getSequence(getNames, as.string=T)), taxaIDs = names1)

## Next, assign taxonomy to genus level with DADA2
set.seed(3874)
taxa <- assignTaxonomy("/Users/mis696/proj/parathaa/input/SILVAsubsample_SeedGenera_V1V2.pcr.fasta", 
                       "/Users/mis696/proj/parathaa/input/20230111.silva.seed_v138_1.ng.dada.fasta",
                       #"/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.seedonly.fa",
                       #"/Users/mis696/proj/parathaa/input/silva_nr99_v138.1_train_set.fa",
                       #"/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.fasta", 
                       multithread=TRUE)
## Remove sequences with undefined ("N") bases, store until after species assignment
taxa.test <- as.data.frame(taxa)
taxa.test$taxaIDs <- names1
nChars <- grep("N", rownames(taxa.test))
print(paste("Removing", length(nChars), "sequences with N bases"))
withNbases <- taxa.test[nChars,]
taxa <- taxa[-nChars,]

## Perform species assignment with DADA2
taxa.sp <- addSpecies(taxa, "/Users/mis696/proj/parathaa/input/silva_species_assignment_v138.1.seedonly.fa")
#                      "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.ng.dada.sp.fasta")

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

ps1_dada.V1V2 <- phyloseq(OTU_dada, TAX_dada, SAMP_dada)

print(ps1_dada.V1V2)


############################################
### Assess performance on Synthetic Data ###
############################################

#Get taxdata from Synthetic.Reads.R
synth.parathaa<- as.data.frame(tax_table(ps1_parathaa.V1V2))
synth.parathaa$AccID <- rownames(synth.parathaa)

synth.parathaa2 <- left_join(synth.parathaa, taxdata, by="AccID")
synth.parathaa2 <- synth.parathaa2 %>% 
  mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
synth.parathaa2 <- synth.parathaa2 %>% 
  mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x))
synth.parathaa3 <- synth.parathaa2 %>% 
  dplyr::rowwise() %>%
  mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
         Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
  )

synth.dada <- as.data.frame(tax_table(ps1_dada.V1V2))
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



##V4V5

test3 <- as.data.frame(tax_table(ps1_parathaa.V4V5))
test3$AccID <- rownames(test3)

test4 <- left_join(test3, taxdata, by="AccID")

test4 <- test4 %>% 
  mutate(Species.x = unlist(lapply(str_split(Species.x, ";"), FUN=function(x) paste0(word(x,1,2), collapse = ";" ))))
test4 <- test4 %>% 
  mutate(Species.x = ifelse(Species.x=="NA", NA, Species.x))
test5 <- test4 %>% 
  dplyr::rowwise() %>%
  mutate(Flag = ifelse(is.na(Species.x), NA, word(Species.y, 1, 2) %in% str_split(Species.x, ";", simplify = T)),
         Flag.genus = ifelse(is.na(Genus.x), NA, Genus.y %in% str_split(Genus.x, ";", simplify = T))
)

test3d <- as.data.frame(tax_table(ps1_dada.V4V5))
test3d$AccID <- rownames(test3d)

test4d <- left_join(test3d, taxdata, by="AccID")
test4d <- test4d %>% mutate(Species.x = word(Species.x, 1, 2))
test5d <- test4d %>% 
  rowwise() %>%
  mutate(Flag = ifelse(is.na(Species.x), 
                       NA, 
                       word(Species.y, 1, 2) %in% t(apply(str_split(Species.x, ";", simplify=T), 1, FUN=function(X) word(X, 1,2)))),
         Flag.genus = ifelse(is.na(Genus.x), 
                             NA, 
                             Genus.y %in% (str_split(Genus.x, ";", simplify=T))
         )
)

compare.synth <- dplyr::full_join(test5d, test5, by="AccID")
compare.synth <- compare.synth %>% 
  mutate(Species.silva = word(Species.y.y, 1, 2)) %>%
  rename(Species.parathaa = Species.x.y,
         Species.dada = Species.x.x,
         Genus.dada = Genus.x.x,
         Genus.parathaa = Genus.x.y,
         Genus.silva = Genus.y.y,
         Phylum.silva = Phylum.y.y)



## Both uniquely correct:
dim(compare.synth %>% filter(Genus.dada==Genus.silva & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
## Paratha uniquely correct, dada2 partly correct:
dim(compare.synth %>% filter(Flag.genus.x & Genus.dada!=Genus.silva & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
## Paratha uniquely correct, dada2 incorrect:
dim(compare.synth %>% filter(!Flag.genus.x & !is.na(Genus.dada)  & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
## Paratha uniquely correct, dada2 unassigned:
dim(compare.synth %>% filter(is.na(Flag.genus.x)  & Genus.parathaa==Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

##Paratha partly correct, dada2 uniquely correct:
dim(compare.synth %>% filter(Genus.dada==Genus.silva & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Both partly correct:
dim(compare.synth %>% filter(Flag.genus.x & Genus.dada!=Genus.silva & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha partly correct, dada2 incorrect:
dim(compare.synth %>% filter(!Flag.genus.x & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha partly correct, dada2 unassigned:
dim(compare.synth %>% filter(is.na(Flag.genus.x) & Flag.genus.y & Genus.parathaa!=Genus.silva) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

##Paratha incorrect, dada2 uniquely correct:
dim(compare.synth %>% filter(Genus.dada==Genus.silva & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha incorrect, dada2 partly correct:
dim(compare.synth %>% filter(Flag.genus.x & Genus.dada!=Genus.silva & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Both incorrect:
dim(compare.synth %>% filter(!Flag.genus.x & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha incorrect, dada2 unassigned:
dim(compare.synth %>% filter(is.na(Flag.genus.x) & !Flag.genus.y ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

##Paratha unassigned, dada2 uniquely correct:
dim(compare.synth %>% filter(Genus.dada==Genus.silva & is.na(Flag.genus.y) ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha unassigned, dada2 partly correct:
dim(compare.synth %>% filter(Flag.genus.x & Genus.dada!=Genus.silva & is.na(Flag.genus.y) ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Paratha unassigned, dada2 incorrect:
dim(compare.synth %>% filter(!Flag.genus.x & is.na(Flag.genus.y)) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))
##Both unassigned:
dim(compare.synth %>% filter(is.na(Flag.genus.x) & is.na(Flag.genus.y) ) %>% select(Genus.dada, Genus.parathaa, Genus.silva, AccID))

###################################
compare.synth <- compare.synth %>% filter(!Species.silva %in% taxdata_SP)
## Both uniquely correct:
dim(compare.synth %>% filter(Species.dada==Species.silva & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
## Paratha uniquely correct, dada2 partly correct:
dim(compare.synth %>% filter(Flag.x & Species.dada!=Species.silva & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
## Paratha uniquely correct, dada2 incorrect:
dim(compare.synth %>% filter(!Flag.x & !is.na(Species.dada)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
## Paratha uniquely correct, dada2 unassigned:
dim(compare.synth %>% filter(is.na(Flag.x)  & Species.parathaa==Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))

##Paratha partly correct, dada2 uniquely correct:
dim(compare.synth %>% filter(Species.dada==Species.silva & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Both partly correct:
dim(compare.synth %>% filter(Flag.x & Species.dada!=Species.silva & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Paratha partly correct, dada2 incorrect:
dim(compare.synth %>% filter(!Flag.x & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Paratha partly correct, dada2 unassigned:
dim(compare.synth %>% filter(is.na(Flag.x) & Flag.y & Species.parathaa!=Species.silva) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))

##Paratha incorrect, dada2 uniquely correct:
dim(compare.synth %>% filter(Species.dada==Species.silva & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Paratha incorrect, dada2 partly correct:
dim(compare.synth %>% filter(Flag.x & Species.dada!=Species.silva & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Both incorrect:
dim(compare.synth %>% filter(!Flag.x & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Paratha incorrect, dada2 unassigned:
dim(compare.synth %>% filter(is.na(Flag.x) & !Flag.y ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))

##Paratha unassigned, dada2 uniquely correct:
dim(compare.synth %>% filter(Species.dada==Species.silva & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Paratha unassigned, dada2 partly correct:
dim(compare.synth %>% filter(Flag.x & Species.dada!=Species.silva & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Paratha unassigned, dada2 incorrect:
dim(compare.synth %>% filter(!Flag.x & is.na(Flag.y)) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))
##Both unassigned:
dim(compare.synth %>% filter(is.na(Flag.x) & is.na(Flag.y) ) %>% select(Species.dada, Species.parathaa, Species.silva, AccID))

## How many assigned 1-to-Many?
length(grep(";", compare.synth$Genus.parathaa))/length(compare.synth$Genus.parathaa)
length(grep(";", compare.synth$Species.parathaa))/length(compare.synth$Species.parathaa)

compare.synth %>% filter(AccID %in% multis$AccID) %>% select(Species.dada, Species.parathaa, Species.y.y, AccID) %>% View
multis %>% filter(Flag.y & Species.parathaa!=Species.silva) %>% View


IncSp <- compare.synth %>% filter(!Flag.y ) %>% select(Species.silva) %>% unique
IncSp <- t(as.matrix(IncSp))
taxdata_seed <- taxdata %>% 
  filter(primaryAccession %in% SeedTax$primaryAccession)
taxdata_SP <- word(taxdata_seed$Species, 1, 2)
## % of incorrect species (unassigned by dada2) that aren't in seed db
1-sum(IncSp[1,] %in% taxdata_SP)/length(IncSp[1,] %in% taxdata_SP)
