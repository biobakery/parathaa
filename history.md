Changes in version 0.2.1 (2025-07-25)
+ Removed nearest neighbor code as the default setting
+ Added in sensitive and specific modes (specific being the default)
+ Added in a new metric for measuring the cloest distance placement when multiple placements occured
+ Added in new default databases using GTDB
+ Added in long read support through the use of "--primers none"
+ Added in a diagnositic script for plotting out sequence placements of assigned sequences
+ Added in support to install via Conda



Changes in version 0.2.1 (2023-11-16)
+ Added nearest neighbor code to improve the specificity of species assignments
+ Began splitting up software pipeline from analysis in paper 
+ Added parrallel support to vastly improve software speed
+ Added a number of fixes to the pipeline to improve speed (mostly using the Castor R package)
+ Added parameters that can tweak over-grouping and over-splitting weights
+ Added parameters to tweak the binomial error model
+ Added support for multiple databases not just SILVA

Changes in version 0.1.1 (2023-04-11):
+ Split make.taxonomy.trees.R into two scripts: find.cutoffs.R and assign.node.tax.R

Changes in version 0.1.0 (2023-03-16):
+ Initial tagged release.
