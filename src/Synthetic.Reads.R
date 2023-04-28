source("src/SILVA.species.editor.R")

## Read in SILVA 138.1 taxonomy for subsetting
inFileTaxdata <- "./input/taxmap_slv_ssu_ref_138.1.txt"

taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
taxdata <- taxdata %>%
  mutate(taxonomy=paste0(path, organism_name))
taxdata <- SILVA.species.editor.DADA(taxdata, "taxonomy")

taxdata <- taxdata %>%
  select(AccID, primaryAccession, start, stop, taxonomy) %>%
  separate(col=taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

## Additional changes (getting rid of subspecies)
taxdata <- SILVA.species.editor(taxdata, Task="find_cutoffs")
taxdata <- taxdata %>% filter(!is.na(Species))

#Use to find species labels with >2 words
#unique(taxdata$Species[sapply(strsplit(taxdata$Species, " ", fixed = TRUE), length)>2])

## Read in seed db for exclusion
inFileSeedDB <- "input/silva.seed_v138_1.tax"
SeedTax <- read.table(inFileSeedDB , header=F, fill=TRUE,sep='\t')

SeedTax <- SeedTax %>%
  separate(col=V1, into=c("primaryAccession", "ArbID"), sep="\\.") %>%
  separate(col=V2, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  filter(Kingdom=="Bacteria" & !is.na(Genus) & Genus!="") 

taxdata_seedless <- taxdata %>% 
  filter(!primaryAccession %in% SeedTax$primaryAccession) %>%
  ## Subset to Genera in Seed for now
  filter(Genus %in% unique(SeedTax$Genus))


## How many to sample from each Genus
GenusCounts <- taxdata_seedless %>% dplyr::count(Genus)
GenusCounts$subsetN <- pmin(pmax(20, ceiling(GenusCounts$n*0.01)), GenusCounts$n)
## Include only genera from seed for now
GenusCounts <- GenusCounts %>% filter(Genus %in% unique(SeedTax$Genus))
dim(GenusCounts)
sum(GenusCounts$subsetN)


library(purrr)
set.seed(978)
subsample <- taxdata_seedless %>% 
  group_split(Genus) %>% 
  map2_dfr(GenusCounts$subsetN, ~ slice_sample(.x, n = .y))

## Print IDs to file
write_lines(subsample$AccID, file="/Users/mis696/proj/parathaa/input/subsampleIDs_SeedGeneraV4V5.txt")

## Subset fasta file with:
## faSomeRecords  /Users/mis696/proj/16s-region-checker/input/SILVA_138.1_SSURef_tax_silva.fasta /Users/mis696/proj/parathaa/input/subsampleIDs.txt /Users/mis696/proj/parathaa/input/SILVAsubsample.fasta 

## faSomeRecords  /Users/mis696/proj/16s-region-checker/input/SILVA_138.1_SSURef_tax_silva.fasta /Users/mis696/proj/parathaa/input/subsampleIDs_SeedGeneraV4V5.txt /Users/mis696/proj/parathaa/input/SILVAsubsample_SeedGenera_V4V5.fasta 

## Trim to specific subregion:
##/Users/mis696/mothur/mothur '#pcr.seqs(fasta = /Users/mis696/proj/parathaa/input/SILVAsubsample_SeedGenera.fasta , oligos=/Users/mis696/proj/parathaa/input/V4V5.oligos.txt, pdiffs=0, rdiffs=0)'

##/Users/mis696/mothur/mothur '#pcr.seqs(fasta = /Users/mis696/proj/parathaa/input/SILVAsubsample_SeedGenera.fasta , oligos=/Users/mis696/proj/parathaa/input/V1V2.oligos.txt, pdiffs=0, rdiffs=0)'

