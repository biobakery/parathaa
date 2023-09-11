#!/usr/bin/env Rscript
require(docopt)
'Usage:
   force_bifurcated_tree.R [-t <tree file> -o <output>]

Options:
   -t tree file in newick format to force bifurcation on
   -o output directory [default: output]

 ]' -> doc

#To add as arguments: binom error params

library(ape)

opts <- docopt(doc)


tree <- read.tree(opts$t)

tree_bifur <- multi2di(tree)

write.tree(tree_bifur, file=file.path(opts$o, "region_specific_bi.tree"))