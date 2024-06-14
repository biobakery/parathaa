#!/usr/bin/env Rscript
require(docopt)
'Usage:
   bind_tax_out.R [-e <exact taxonomy file> -t <normal taxonomy file>]

Options:
   -e exact taxonomies
   -t taxonomies

 ]' -> doc
opts <- docopt(doc)


exact_taxonomy <- read.table(opts$e, sep="\t", header=T)
norm_taxonomy <- read.table(opts$t, sep="\t", header=T)


merged_taxonomy <- rbind(exact_taxonomy[,-10], norm_taxonomy)

write.table(merged_taxonomy, file=opts$t, sep="\t", quote=F, row.names=F)
