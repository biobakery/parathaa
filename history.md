Changes in version 0.1.0 (2023-03-16):
+ Initial tagged release.

Changes in version 0.1.1 (2023-04-11):
+ Split make.taxonomy.trees.R into two scripts: find.cutoffs.R and assign.node.tax.R

Changes in version 0.2.1 (2023-11-16)
+ Added nearest neighbor code to improve the specificity of species assignments
+ Began splitting up software pipeline from analysis in paper 
+ Added parrallel support to vastly improve software speed
+ Added a number of fixes to the pipeline to improve speed (mostly using the Castor R package)
+ Added parameters that can tweak over-grouping and over-splitting weights
+ Added parameters to tweak the binomial error model
+ Added support for multiple databases not just SILVA
