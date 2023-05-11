# convert Kelsey's data to fasta
kelseyData <- read.delim("~/Downloads/all_samples_taxonomy_closed_reference_withseqs.tsv", sep='\t', fill=T, stringsAsFactors = F)

kelseyData <- kelseyData[c(1,126)] 
taxaIDs <- paste("ID", 1:nrow(kelseyData), sep = "")
kelseyData <- cbind(taxaIDs, kelseyData[,2], kelseyData[,1])
kelseyData <- as.data.frame(kelseyData)
colnames(kelseyData) <- c("taxaIDs", "Taxonomy_1", "seq")

#View(kelseyData)


Xfasta <- character(nrow(kelseyData) * 2)
Xfasta[c(TRUE, FALSE)] <- paste0(">", kelseyData[,"taxaIDs"], '\t', kelseyData[,"Taxonomy_1"])
Xfasta[c(FALSE, TRUE)] <- kelseyData[,"seq"]
#View(Xfasta)

writeLines(Xfasta, "~/proj/PPITAA/input/amplicons_fromKelsey.fasta")




