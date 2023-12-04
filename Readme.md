# ***ATTENTION***

Before opening a new issue here, please check the appropriate help channel on the bioBakery Support Forum (https://forum.biobakery.org) and consider opening or commenting on a thread there.

----

# PARATHAA User Manual 
> Preserving and Assimilating Region-specific Ambiguities in Taxonomic Hierarchical Assignments for Amplicons

----


**If you use PARATHAA in your work, please cite PARATHAA at:**
[http://huttenhower.sph.harvard.edu/PARATHAA](http://huttenhower.sph.harvard.edu/PARATHAA)

----


**For additional information, read the [PARATHAA Tutorial](https://github.com/biobakery/biobakery/wiki/PARATHAA)**

PARATHAA(Preserving and Assimilating Region-specific Ambiguities in Taxonomic Hierarchical Assignments for Amplicons) is a tool used for the taxonomic assignment of 16S rRNA gene sequences that takes into account the uncertainty associated with using specific variable regions/primers. PARATHAA does this by generating new primer-trimmed phylogenetic trees from reference 16S rRNA gene datasets and then determines the optimal phylogenetic distances within that tree for taxonomic labeling. PARATHAA then can use this tree to assign taxonomy to query 16S rRNA gene sequences by aligning and placing those sequences into the new primer-trimmed reference database.

----

## Contents 
* [Features](#Features)
* [Requirements](#Requirements)
* [Installation](#Installation)
  * [Manual Installation](#Manual_Installation)
  * [Pypi Installation](#Pypi_Installation)
* [Workflow](#workflow)
  * [Step1](#Step1)
  * [Step2](#Step2)   
* [Databases](#Databases)
* [Outputs](#Outputs)
* [Complete option list](#Complete_option_list)


----


####  Manual Installation
- Clone your recently created repository in your local development environment either using:
```
git clone https://github.com/biobakery/parathaa
```
or using the "**Clone or Download**" button.

#### Requirements
To run Parathaa, you should be running linux/MacOS, and you should have installed:
- Python >=3.7
- R > 4.0.0
 - With the following R libraries
 - Castor (https://cran.r-project.org/web/packages/castor/index.html)
 - Stringr (https://cran.r-project.org/web/packages/stringr/index.html)
 - docopt (https://cran.r-project.org/web/packages/docopt/index.html)
 - ggtree (https://bioconductor.org/packages/release/bioc/html/ggtree.html)
 - tidytree (https://cran.r-project.org/web/packages/tidytree/index.html)
 - treeio (https://bioconductor.org/packages/release/bioc/html/treeio.html)
 - dplyr (https://cran.r-project.org/web/packages/dplyr/index.html)
 - phytools (https://cran.r-project.org/web/packages/phytools/index.html)
 - doSNOW (https://cran.r-project.org/web/packages/doSNOW/index.html)
 - tidyr (https://cran.r-project.org/web/packages/tidyr/index.html)
 - ape (https://cran.r-project.org/web/packages/ape/index.html)
 - TDbook (https://cran.r-project.org/web/packages/TDbook/index.html)
 - readr (https://cran.r-project.org/web/packages/readr/index.html)
 - ggplot2 (https://cran.r-project.org/web/packages/ggplot2/index.html)
 - reshape2 (https://cran.r-project.org/web/packages/reshape2/index.html)
- AnADAMA2: https://pypi.org/project/anadama2/
- mothur: https://anaconda.org/bioconda/mothur
- taxtastic: https://pypi.org/project/taxtastic/
- fasttree: http://www.microbesonline.org/fasttree/#Install
- pplacer: https://anaconda.org/bioconda/pplacer


#### Workflow:
- Parathaa is seperated into two seperate workflows that are used in conjunction with one another. In general most users will only use step 2 with already pre-computed files from step 1. However users able to adjust the reference database parathaa uses for taxonomic assignment using the commands in step 1.

#### Step 1: Creation of primer trimmed phylogenetic trees with taxonomic labels

- The first step in parathaa is to create a primer trimmed phyloegentic tree that has its internal nodes labelled with taxonomic annotations. This is done using the command:
```
run_tree_analysis.py
```

In brief this command does the following:

1. Takes in the reference MSA and trims the sequences to the region that was amplified based on the given set of primers
2. Generates a new phylogenetic tree based on these newly trimmed sequences using FastTree
3. Finds the appropriate distance thresholding cutoffs for each taxonomic level
4. Assigns taxonomy to the internal nodes of the new primer trimmed reference tree


#### Inputs
This command takes in the following:
- `primers` A set of primers.
- `database` A reference database that consists of a pre-aligned 16S rRNA gene reference set. By default we a use silva v138 seed database.
- `taxonomy` A taxonomy file that contains the original taxonomic labels for the 16S sequences in the provided reference database

There are two optional commands as well:
- ```swight``` which controls the penalty weighting for over splitting taxonomic groups
- ```mwight``` which controls the penalty weighting for over merge taxonomic groups

These two options are used when calculating the optimal threshold distance along the tree for each taxonomic level. 

##### Demo Run Step 1:

Not some input files will have been gzip'd so that they can fit within this repository. Please check the input files before running and expand them with ```gunzip```. 
```
run_tree_analysis.py --primers input/primers/V4V5.oligos --database input/testing/tree_construction/subset_silva_seed_v138_1.align  --output output --taxonomy input/silva_v138/taxmap_slv_ssu_ref_138.1.txt
```

#### Step 2: Taxonomic assignment

PARATHAA will come with a number of pre-computed-trees for this step so that users do not need to generate their own primer trimmed trees for commonly used reference databases. However this is still under development and currently the only pre-computed-tree available is for silva_v138 using V4V5 primers. 


- The second step of PARATHAA uses the phylogenetic tree, multiple sequence alignment (MSA) and tree data files to create taxonomic assignments to 16S rRNA gene sequences. This step is running using the following command:
```
run_taxa_assignment.py
```

Briefly this step takes in the newly created primer trimmed tree, MSA, tree reference files, and thresholding information to assign taxonomy to query 16S sequences by the following steps.

1. Aligns query sequences to primer trimmed MSA generated in step 1
2. Places those sequences into the primer trimmed tree generated in step 1 using pplacer
3. Assigns taxonomy to those query sequences based on their placement and their distance from taxonomically labelled interior nodes. Note the thresholding distance computed in step one determines the apprioriate distance away a node can be to be assigned to a taxonomic level.

##### Inputs:
This command takes the following inputs:
- `trimmedDatabase` the primer trimmed MSA generated from step 1
- `trimmedTree` the trimmed phylogenetic tree generated from step 1
- `treeLog` the treelog.txt file generated from step 1
- `query` the 16S sequences that you want to assign taxonomy
- `thresholds` an RData file containing the optimal phylogenetic distances for taxonomic assignment
- `namedTree` an RData file containing the annotated trimmed phylogenetic tree  

##### Demo Run Step 2:
Not some input files will have been gzip'd so that they can fit within this repository. Please check the input files before running and expand them with ```gunzip```.

```
run_taxa_assignment.py --trimmedDatabase input/silva_v138/Pre-computed-Trees/V4V5/silva.seed_v138_1.pcr.align \
--trimmedTree input/silva_v138/Pre-computed-Trees/V4V5/region_specific.tree \
--treeLog input/silva_v138/Pre-computed-Trees/V4V5/treelog.txt \
--query input/testing/taxa_assignment/SRR3225703_V4V5_subset.fasta \
--thresholds input/silva_v138/Pre-computed-Trees/V4V5/optimal_scores.RData \
--namedTree input/silva_v138/Pre-computed-Trees/V4V5/resultTree_bestThresholds.RData \
--output output_taxa_test
```

