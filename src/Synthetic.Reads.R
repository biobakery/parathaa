## Read in SILVA NR99 taxonomy for subsetting
inFileTaxdata <- "/Users/mis696/proj/16s-region-checker/input/taxmap_slv_ssu_ref_138.1.txt"

taxdata <- read.table(inFileTaxdata , header=T, fill=TRUE,sep='\t', quote="")
taxdata <- taxdata %>%
  unite("AccID", c("primaryAccession", "start", "stop"), sep=".", remove=F)
taxdata <- taxdata %>%
  separate(col=path, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  dplyr::rename(Species = organism_name) %>%
  filter(Kingdom=="Bacteria") 
taxdata <- taxdata %>% 
  filter(!is.na(Genus) & Genus!="")

## Read in seed db for exclusion
inFileSeedDB <- "/Users/mis696/proj/parathaa/input/silva.seed_v138_1.tax"
SeedTax <- read.table(inFileSeedDB , header=F, fill=TRUE,sep='\t')
SeedTax <- SeedTax %>%
  separate(col=V1, into=c("primaryAccession", "ArbID"), sep="\\.") %>%
  separate(col=V2, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep=";") %>%
  filter(Kingdom=="Bacteria" & !is.na(Genus) & Genus!="") 

taxdata_seedless <- taxdata %>% 
  filter(!primaryAccession %in% SeedTax$primaryAccession)
             
## How many to sample from each Genus
GenusCounts <- taxdata_seedless %>% dplyr::count(Genus)
GenusCounts$subsetN <- pmin(pmax(5, ceiling(GenusCounts$n*0.01)), GenusCounts$n)
sum(GenusCounts$subsetN)

library(purrr)

subsample <- taxdata_seedless %>% 
  group_split(Genus) %>% 
  map2_dfr(GenusCounts$subsetN, ~ slice_sample(.x, n = .y))

## Print IDs to file
write_lines(subsample$AccID, file="/Users/mis696/proj/parathaa/input/subsampleIDs.txt")

## Subset fasta file with:

## faSomeRecords  /Users/mis696/proj/16s-region-checker/input/SILVA_138.1_SSURef_tax_silva.fasta /Users/mis696/proj/parathaa/input/subsampleIDs.txt /Users/mis696/proj/parathaa/input/SILVAsubsample.fasta 

## Trim to specific subregion:
## /Users/mis696/mothur/mothur '#pcr.seqs(fasta = /Users/mis696/proj/parathaa/input/SILVAsubsample.fasta , oligos=/Users/mis696/proj/16s-region-checker/input/EMPV4.oligos, pdiffs=0, rdiffs=0)'


