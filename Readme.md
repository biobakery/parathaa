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

PARATHAA (Preserving and Assimilating Region-specific Ambiguities in Taxonomic Hierarchical Assignments for Amplicons) is a tool used for the taxonomic assignment of 16S rRNA gene sequences that takes into account the uncertainty associated with using specific variable regions/primers. PARATHAA does this by generating new primer-trimmed phylogenetic trees from reference 16S rRNA gene datasets and then determines the optimal phylogenetic distances within that tree for taxonomic labeling. PARATHAA then can use this tree to assign taxonomy to query 16S rRNA gene sequences by aligning and placing those sequences into the new primer-trimmed reference database.

----

## Contents 
* [Features](#Features)
* [Requirements](#Requirements)
* [Installation](#Installation)
  * [Manual Installation](#Manual_Installation)
  * [Pypi Installation](#Pypi_Installation)
* [Databases](#Databases)
* [Workflow](#workflow)
  * [Step1](#Step1)
  * [Step2](#Step2)   
* [Diagnostics](#PARATHAA_Diagnostic_Scripts)


----
### Requirements

Below are the required dependencies to run Parathaa, you should be running linux/macOS (note if running macOS follow special instructions for pplacer install). We have included installation commands for various dependencies  If the user wants to install outside of conda please check out each tool's main page for other installation instructions.

For our next release we plan to create a conda package that will cover all of these dependencies automatically.

- Conda (we suggest miniconda https://docs.conda.io/projects/miniconda/en/latest/ although any flavor should work)
- Python >=3.7
- R > 4.0.0

#### Create conda environment and set up R and python

Create a new conda environment and activate it:

```
conda create -n parathaa.0.1 python=3.9

conda activate parathaa.0.1
```

Setup conda with R

```
conda install -c conda-forge r-base=4.2.3
conda install -c r r-openssl
```

Load R, install packages, quit:

```
R
```

```

options(timeout=600)
install.packages(c("ape", 
                   "castor",
                   "stringr", 
                   "docopt",
                   "tidytree",
                   "dplyr",
                   "phytools",
                   "doSNOW",
                   "tidyr",
                   "TDbook",
                   "ggplot2",
                   "reshape2")) # add Ncpus = 4 to go faster

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("ggtree",
                       "treeio"))

quit()
```

- AnADAMA2: https://pypi.org/project/anadama2/

```
pip install anadama2
```

- mothur: https://anaconda.org/bioconda/mothur

```
conda install -c bioconda mothur
```

- taxtastic: https://pypi.org/project/taxtastic/

```
pip install taxtastic
```

- fasttree: http://www.microbesonline.org/fasttree/#Install

```
conda install -c bioconda fasttree
```

Or install via their website. We recommend using the multithreaded executable which is supported by parathaa. After downloading the executable please add it to PATH within your computing environment.
 
- pplacer: https://anaconda.org/bioconda/pplacer

**If using linux:**

```
conda install -c bioconda pplacer
```

**If using macOS please install SEPP which includes a recent compiled version of pplacer for macOSX. (we will package this up into its own conda environemnt in the future)**

```
conda install -c bioconda sepp
```


###  Manual Installation of PARATHAA
- Clone your recently created repository in your local development environment either using:
```
git clone https://github.com/biobakery/parathaa
```
or using the "**Clone or Download**" button.

- install scripts to your bin directory of choice

```
python setup.py install 
```

### Installation of PARATHAA through pip

```
pip install parathaa 
```
**NOTE:**
If you do not have write permissions to '/usr/lib/', then add the option "--prefix=<CUSTOM_PYTHON_LIBS_PATH>" to the parathaa install command. This will install the python package into subdirectories of '<CUSTOM_PYTHON_LIBS_PATH>' on Linux. Please note when using the "--prefix" install option on some platforms, you might need to add 'CUSTOM_PYTHON_LIBS_PATH' to your $PYTHONPATH as it might not be included by default. 
```
export PYTHONPATH=$PYTHONPATH:<CUSTOM_PYTHON_LIBS_PATH>
```


### Databases

As part of PARATHAA we provide pre-computed databases for commonly used 16S variable regions on a SEED database of SILVA (provided by the mothur team here: https://mothur.org/wiki/silva_reference_files/)

We are currently working with version 138.1

Use the following command to download the Pre-computed databases:


16S Full Length - SILVA SEED 138.1

```
wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_FL.tar.gz
tar -xvf SILVA_FL.tar.gz
rm SILVA_FL.tar.gz
```

16S V1V2 - SILVA SEED 138.1

```
wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_V1V2.tar.gz
tar -xvf SILVA_V1V2.tar.gz
rm SILVA_V1V2.tar.gz
```

16S V3V4 - SILVA SEED 138.1

```
wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_V3V4.tar.gz
tar -xvf SILVA_V3V4.tar.gz
rm SILVA_V3V4.tar.gz
```

16S V4V5 - SILVA SEED 138.1

```
wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_V4V5.tar.gz
tar -xvf SILVA_V4V5.tar.gz
rm SILVA_V4V5.tar.gz
```


### Workflow:
- Parathaa is seperated into two seperate workflows that are used in conjunction with one another. In general most users will only use step 2 with already pre-computed files from step 1. However users able to adjust the reference database parathaa uses for taxonomic assignment using the commands in step 1.

#### Step 1: Creation of primer trimmed phylogenetic trees with taxonomic labels

- The first step in parathaa is to create a primer trimmed phyloegentic tree that has its internal nodes labelled with taxonomic annotations. This is done using the command:
```
parathaa_run_tree_analysis
```

In brief this command does the following:

1. Takes in the reference MSA and trims the sequences to the region that was amplified based on the given set of primers
2. Generates a new phylogenetic tree based on these newly trimmed sequences using FastTree
3. Finds the appropriate distance thresholding cutoffs for each taxonomic level
4. Assigns taxonomy to the internal nodes of the new primer-trimmed reference tree


#### Inputs
This command takes in the following:
- `primers` A set of primers.
- `database` A reference database that consists of a pre-aligned 16S rRNA gene reference set. By default we a use silva v138 seed database.
- `taxonomy` A taxonomy file that contains the original taxonomic labels for the 16S sequences in the provided reference database

There are two optional commands as well:
- ```sweight``` which controls the penalty weighting for over splitting taxonomic groups
- ```mweight``` which controls the penalty weighting for over merge taxonomic groups

These two options are used when calculating the optimal threshold distance along the tree for each taxonomic level, and their default values are each 1, which weights over-splitting and over-merging equally. 

##### Demo Run Step 1:

Not some input files will have been gzip'd so that they can fit within this repository. Please check the input files before running and expand them with ```gunzip```. 

```
gunzip input/silva_v138/taxmap_slv_ssu_ref_138.1.txt.gz
```

```
parathaa_run_tree_analysis --primers input/primers/V4V5.oligos --database input/testing/tree_construction/subset_silva_seed_v138_1.align  --output output --taxonomy input/silva_v138/taxmap_slv_ssu_ref_138.1.txt
```

#### Step 2: Taxonomic assignment

PARATHAA will come with a number of pre-computed-trees for this step so that users do not need to generate their own primer trimmed trees for commonly used reference databases. However this is still under development and currently the only pre-computed-tree available is for silva_v138 using V4V5 primers. 


- The second step of PARATHAA uses the phylogenetic tree, multiple sequence alignment (MSA) and tree data files to create taxonomic assignments to 16S rRNA gene sequences. This step is running using the following command:
```
parathaa_run_taxa_assignment
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

There are two optional commands as well, related to species-level classifications:
- ```delta``` which controls which species assignments are removed due to nearby reference sequences from different species. Default value is 0.005; we recommend setting this to 0.5 * the species threshold identified by Parathaa for your region in step 1.
- ```mult``` which shrinks the species threshold for better performance. Default value is 0.5, which is recommended.

##### Demo Run Step 2:


For this demo we will first need to download the Pre-computed V4V5 database.

```
wget  http://huttenhower.sph.harvard.edu/parathaa_db/SILVA_V4V5.tar.gz
tar -xvf SILVA_V4V5.tar.gz
rm SILVA_V4V5.tar.gz
```


```
parathaa_run_taxa_assignment --trimmedDatabase SILVA_V4V5/silva.seed_v138_1.pcr.align --treeLog SILVA_V4V5/treelog.txt --query /input/testing/taxa_assignment/SRR3225703_V4V5_subset.fasta --thresholds SILVA_V4V5/optimal_scores.RData --namedTree SILVA_V4V5/resultTree_bestThresholds.RData --output test --trimmedTree SILVA_V4V5/region_specific.tree
```

OR

Alternatively you can give this command a directory that contains the ```trimmedTree```, ```treeLog```, ```trimmedDatabase```, ```thresholds```, ```namedTree```, files using the parameter ```--treeFiles```.

```
parathaa_run_taxa_assignment --treeFiles SILVA_V4V5/ \
--query input/testing/taxa_assignment/SRR3225703_V4V5_subset.fasta \
--output output_taxa_test
```

### PARATHAA Diagnostic Scripts


- plot_placement.R

This script contains a useful function for plotting placements. The function is called by other diagnostic scripts. Its important that the input tree has unique labels for all tips and nodes. This can be accomplished with the following R code:

```
in.tree$label <- make.unique(in.tree$label)
```


