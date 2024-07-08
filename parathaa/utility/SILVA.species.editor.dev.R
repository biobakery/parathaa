## Manual/automatic SILVA species consolidation
library(dplyr)
library(stringr)
library(tidyr)

SILVA.species.editor <- function(dataf){

spNames <- unique(dataf$Species)
spNames <- str_replace_all(spNames, pattern = "\\[", replacement = "")
spNames <- str_replace_all(spNames, pattern = "\\]", replacement = "") 

## remove if "bacterium" at start <- strain-specific
spNames2.1 <- spNames[startsWith(spNames, "bacterium")]

spNamesdf <- data.frame(spNames) 
spNamesdf <- spNamesdf %>% mutate(spNames2.1 =replace(spNames, 
                                       spNames %in% spNames2.1, 
                                       NA))

## remove "subsp" <- strain-specific within species? i.e. keep before
spNames2 <-  grep("subsp\\.", spNamesdf$spNames2.1, value=T)
spNames2replace <- str_split(spNames2, " subsp\\.",simplify = T)[,1]

spNamesdf <- spNamesdf %>% mutate(spNames2 =replace(spNames2.1, 
                    spNames %in% spNames2, 
                    spNames2replace))

## remove if "bacterium" <- mix of specific strains and general "uncultured XXX bacterium"
spNames3 <- grep(" bacterium", spNamesdf$spNames2, value = T)

spNamesdf <- spNamesdf %>% mutate(spNames3 =replace(spNames2, 
                                                    spNames %in% spNames3, 
                                                    NA))

## remove if "uncultured" <- mix of specific strains and general "uncultured proteobacterium"
spNames4 <- grep("uncultured", spNamesdf$spNames3, value = T)

spNamesdf <- spNamesdf %>% mutate(spNames4 =replace(spNames3, 
                                                    spNames %in% spNames4, 
                                                    NA))

## remove if "metagenome" <- general, only 1
spNames4.1 <- grep("metagenome", spNamesdf$spNames4, value = T)

spNamesdf <- spNamesdf %>% mutate(spNames4.1 =replace(spNames4, 
                                                    spNames %in% spNames4.1, 
                                                    NA))
# remove if contains "sp." <- mix of general and specific. ends with sp. are general, can try to remove those first, separately
spNames5 <- grep("sp\\.", spNamesdf$spNames4.1, value = T)

spNamesdf <- spNamesdf %>% mutate(spNames5 =replace(spNames4.1, 
                                                    spNames %in% spNames5, 
                                                    NA))


# remove if "unidentified" <- mix of general ("unidentified", "unidentified proteobacterium") and specific ("unidentified proteobacterium v1")
spNames6 <- grep("unidentified", spNamesdf$spNames5, value = T)

spNamesdf <- spNamesdf %>% mutate(spNames6 =replace(spNames5, 
                                                    spNames %in% spNames6, 
                                                    NA))

# remove unidentified strains <- actually seem to be identified?
spNames7 <- spNamesdf$spNames6[str_split(spNamesdf$spNames6, " ", simplify = T)[,2]=="str."]
spNamesdf <- spNamesdf %>% mutate(spNames7 =replace(spNames6, 
                                                    spNames %in% spNames7, 
                                                    NA))

# recode identified strains
spNames8 <- grep("str\\.", spNamesdf$spNames7, value = T)
spNames8replace <- str_split(spNames8, " ",simplify = T)[,1:2]
spNames8replace <- paste(spNames8replace[,1], spNames8replace[,2])

spNamesdf <- spNamesdf %>% mutate(spNames8 =replace(spNames7, 
                                                    spNames %in% spNames8, 
                                                    spNames8replace))

nameTest <- spNames[grep(pattern = "^[A-Z][a-z]+ [a-z]+$" ,x = spNames, perl=TRUE, value = F)]
spNames[!spNames %in% nameTest]


## Code out strain designations for species
spNames9 <- spNamesdf$spNames8[str_split(spNamesdf$spNames8, " ",simplify = T)[,3]!=""]

Index <- c("[Haemophilus] ducreyi", "Acidiphilium cryptum", "Acidothermus cellulolyticus", "Actinobacillus capsulatus",
           "Actinobacillus pleuropneumoniae", "Afipia felis", "Agrobacterium radiobacter", 
           "Agrobacterium vitis", "Aliivibrio fischeri", "Alkalilimnicola ehrlichii",
           "Aphanizomenon flos-aquae", "Aquifex aeolicus", "Bacillus amyloliquefaciens", 
           "Bacillus cereus", "Bacillus halodurans", "Bacillus lehensis",
           "Bacillus subtilis", "Bacillus thuringiensis", "Bacillus velezensis",
           "Bacteroides fragilis", "Borreliella burgdorferi", "Bordetella bronchiseptica", "Bradyrhizobium diazoefficiens", "Brucella ceti", "Brucella suis", 
           "Buchnera aphidicola", "Butyrivibrio crossotus", "Caldora penicillata",
           "Calothrix brevissima", "Capnocytophaga sputigena", "Caulobacter segnis", "Caulobacter vibrioides",
           "Cellulomonas flavigena", "Cellulomonas massiliensis", "Cellulophaga algicola",
           "Cephalothrix komarekiana", "Chitinophaga pinensis", "Chlamydia caviae", 
           "Chlamydia trachomatis", "Chlorobium phaeobacteroides", "Chryseobacterium hispalense",
           "Clostridium acetobutylicum", "Clostridium botulinum", "Clostridium perfringens",
           "Clostridium senegalense", "Corynebacterium glutamicum", "Croceibacter atlanticus",
           "Cylindrospermum stagnale", "Desulfococcus oleovorans", "Desulfotalea psychrophila",
           "Desulfovibrio vulgaris",  "Dichelobacter nodosus", "Dickeya chrysanthemi", 
           "Dolichospermum planctonicum", "Edwardsiella ictaluri",
           "Enterococcus faecalis", "Escherichia coli", "Fischerella major",
           "Fluviicola taffensis", "Geminocystis herdmanii", "Geobacter sulfurreducens",
           "Gordonia alkanivorans", "Haemophilus influenzae", "Hahella chejuensis",
           "Halothiobacillus neapolitanus", "Helicobacter hepaticus", "Helicobacter pylori",
           "Henriciella marina", "Histophilus somni", "Hydrocarboniphaga effusa",
           "Idiomarina loihiensis", "Komagataeibacter xylinus", "Lactobacillus casei", 
           "Lactobacillus crispatus",  "Lactobacillus gasseri", "Lactobacillus johnsonii",
           "Lactobacillus plantarum", "Lactobacillus salivarius", "Lactococcus garvieae",
           "Legionella longbeachae", "Limnothrix redekei", "Listeria monocytogenes",
           "Loriellopsis cavernicola", "Marinobacter hydrocarbonoclasticus", "Mesorhizobium japonicum",
           "Mesorhizobium jarvisii", "Methylobacterium nodulans", "Methylocella silvestris",
           "Methyloglobulus morosus", "Methylorubrum extorquens", "Microcystis aeruginosa", "Microcystis wesenbergii",
           "Muricauda ruestringensis", "Mycobacterium tuberculosis", "Mycolicibacterium smegmatis",
           "Myxococcus fulvus",  "Myxococcus xanthus", "Neisseria meningitidis", "Nitrosococcus oceani", "Nitrospira marina",
           "Nostocoida limicola", "Novosphingobium aromaticivorans", "Oceanobacillus iheyensis",
           "Ochrobactrum intermedium", "Oscillatoria nigro-viridis", "Paraglaciecola mesophila",
           "Pasteurella aerogenes", "Pasteurella haemolytica", "Pectobacterium atrosepticum", "Pelobacter propionicus",
           "Phormidium ambiguum", "Polaribacter irgensii", "Pseudanabaena catenata",
           "Pseudoalteromonas atlantica", "Pseudomonas aeruginosa", "Pseudomonas amygdali",
           "Pseudomonas coronafaciens", "Pseudomonas protegens", "Pseudomonas putida",
           "Pseudomonas savastanoi", "Pseudomonas stutzeri", "Pseudomonas syringae",
           "Pseudooceanicola batsensis", "Pseudopropionibacterium propionicum", "Psychromonas ingrahamii", "Rhizobium leguminosarum",
           "Rhodobacter sphaeroides", "Rickettsia bellii", "Rubrobacter xylanophilus",
           "Salinispora arenicola", "Shewanella oneidensis", "Shewanella pealeana", 
           "Shigella boydii", "Shigella flexneri", "Sinorhizobium fredii", "Sinorhizobium meliloti",
           "Sorangium cellulosum", "Sphingopyxis alaskensis", "Staphylococcus aureus",
           "Staphylococcus epidermidis", "Streptococcus agalactiae", "Streptococcus mutans",
           "Streptococcus pneumoniae",
           "Streptococcus salivarius", "Streptomyces ambofaciens", "Streptomyces avermitilis",
           "Streptomyces bottropensis", "Streptomyces fulvissimus", "Streptomyces globisporus",
           "Streptomyces lividans", "Streptomyces scabiei", "Sulfurospirillum deleyianum", 
           "Symploca atlantica", "Synechococcus elongatus", "Syntrophus aciditrophicus",
           "Thermanaerovibrio acidaminovorans", "Thermosynechococcus elongatus", "Thermus thermophilus",
           "Tolumonas auensis", "Trichormus variabilis", "Vibrio cholerae", "Vibrio vulnificus", "Xanthomonas campestris",
           "Xanthomonas citri", "Xanthomonas phaseoli", "Xylella fastidiosa", 
           "Yersinia pestis", "Yersinia pseudotuberculosis")
           
spNamesdf$spNames9 <- spNamesdf$spNames8
for(ind in Index){
spNamesdf <- spNamesdf %>% mutate(spNames9 =replace(spNames9, 
                                                    grep(ind, spNames), 
                                                    ind))
}

## Recode species names in dataset

renamingVec <- spNamesdf$spNames9
names(renamingVec) <- spNamesdf$spNames

newNames <- renamingVec[dataf$Species]

dataf$Species <- newNames

## remove cases where genus name doesn't match Species
spNames10 <- strsplit(dataf$Species, split= " ")

spGenera <- sapply(spNames10, function(x) x[1])

remove <- !str_detect(dataf$Genus,spGenera)

dataf$Species[remove] <- NA

#remove uncultured higher-level taxa
unc.class <- str_detect(dataf$Class, "uncultured") 
dataf$Class[unc.class]<-NA
dataf$Order[unc.class]<-NA
dataf$Family[unc.class]<-NA
dataf$Genus[unc.class]<-NA
dataf$Species[unc.class]<-NA

unc.order <- str_detect(dataf$Order, "uncultured") 
dataf$Order[unc.order]<-NA
dataf$Family[unc.order]<-NA
dataf$Genus[unc.order]<-NA
dataf$Species[unc.order]<-NA

unc.family <- str_detect(dataf$Family, "uncultured") 
dataf$Family[unc.family]<-NA
dataf$Genus[unc.family]<-NA
dataf$Species[unc.family]<-NA

unc.genus <- str_detect(dataf$Genus, "uncultured") 
dataf$Genus[unc.genus]<-NA
dataf$Species[unc.genus]<-NA

return(dataf)
}

## matchGenera function taken from DADA2
matchGenera <- function(gen.tax, gen.binom, split.glyph="/") {
  if(is.na(gen.tax) || is.na(gen.binom)) { return(FALSE) }
  if((gen.tax==gen.binom) || 
     grepl(paste0("^", gen.binom, "[ _", split.glyph, "]"), gen.tax) || 
     grepl(paste0(split.glyph, gen.binom, "$"), gen.tax)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#dataf <- dataf %>% mutate(taxfordada = paste0(path, organism_name))
#colname <- "taxfordada"

SILVA.species.editor.DADA <- function(dataf, colname){
  is.bact <- grepl("Bacteria;", dataf[,colname], fixed=TRUE)
  dataf <- dataf[is.bact,]
  is.uncult <- grepl("[Uu]ncultured", dataf[,colname])
  dataf <- dataf[!is.uncult,]
  is.unident <- grepl("[Uu]nidentified", dataf[,colname])
  dataf <- dataf[!is.unident,]
  is.complete <- sapply(strsplit(as.character(dataf[,colname]), ";"), length)==7
  dataf <- dataf[is.complete,]
  
  # Pull out binomial strings
  binom <- strsplit(as.character(dataf[,colname]), ";")
  genus <- sapply(binom, `[`, 6)
  binom <- sapply(binom, `[`, 7)
  
  genus <- gsub("Candidatus ", "", genus)
  binom <- gsub("Candidatus ", "", binom)
  genus <- gsub("\\[", "", genus)
  genus <- gsub("\\]", "", genus)
  binom <- gsub("\\[", "", binom)
  binom <- gsub("\\]", "", binom)
  
  # Subset down to those binomials which match the curated genus
  genus.binom <- sapply(strsplit(binom, "\\s"), `[`, 1)
  gen.match <- mapply(matchGenera, genus, genus.binom, split.glyph="-")
  # Note that raw Silva files use Escherichia-Shigella, but this is changed to Escherichia/Shigella in dada2 version
  dataf <- dataf[gen.match,]
  binom <- binom[gen.match]
  genus <- genus[gen.match]
  
  # Make matrix of genus/species
  binom[sapply(strsplit(binom, "\\s"), length)==1] <- paste(binom[sapply(strsplit(binom, "\\s"), length)==1], "sp.")
  binom2 <- cbind(sapply(strsplit(binom, "\\s"), `[`, 1),
                  sapply(strsplit(binom, "\\s"), `[`, 2))
  # Keep only those with a named species
  has.spec <- !grepl("sp\\.", binom2[,2]) & !(binom2[,2]=="endosymbiont")
  binom2 <- binom2[has.spec,]
  dataf <- dataf[has.spec,]
  binom <- binom[has.spec]
  genus <- genus[has.spec]
  cat(length(binom), "sequences with genus/species binomial annotation output.\n")
  return(dataf)
}