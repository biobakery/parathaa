#!/usr/bin/python3

### This script given an aligned 16S file, an accompanying taxonomy for those sequences
### and a primer pair will output a phylogenetic tree based on the regions amplified by those primers


#### TO DO:

# allow users to to skip primer trimming step if they just want full length and only the thresholding...
# set up the output so that they go into their own folder...
# generate a small test case...


import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable
from pathlib import Path


# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.2.1",
    description=("Preserving Primer-Induced Taxonomic Ambiguities for Amplicons: This script given an aligned 16S file, an accompanying taxonomy file "
                 " and a primer pair will output a phylogenetic tree of those 16S sequences based on the amplified regions and find optimal "
                 " threhold distances for taxonomy assignment and assigns taxonomy to the internal nodes of the tree"   #Update the description as needed
    ) 
)

# Setting additional custom arguments for workflow - run_tree_analysis
workflow.add_argument(
    name="primers",
    desc="File with primer oligos [default: input/EMPV4.oligos]",
    default="input/EMPV4.oligos")
    
workflow.add_argument(
    name="database",
    desc="Database for taxonomy which is represented by aligned 16S sequences [default: input/silva.seed_v138_1/silva.seed_v138_1.align]",
    default="input/silva.seed_v138_1/silva.seed_v138_1.align")

workflow.add_argument(
    name="taxonomy",
    desc="File that contains taxonomy for the input database",
    default="input/taxmap_slv_ssu_ref_138.1.txt"
)

workflow.add_argument(
    name="sweight",
    desc="penalty weight for over-splitting errors [default: 1]",
    default=1
)

workflow.add_argument(
    name="mweight",
    desc="penalty weight for over-merging errors [default: 1]",
    default=1
)

workflow.add_argument(
    name="errorRate",
    desc="The error rate parameter for the binomial model for ambigious taxonomic assignment",
    default=0.05
)

workflow.add_argument(
    name="binoThreshold",
    desc="Threshold for the binomial error model for ambigious taxonomic assignment",
    default=0.20
)


workflow.add_argument(
    name="threads",
    desc="Number of threads to run multi-threaded processes",
    default="1"

)

# Parsing the workflow arguments
args = workflow.parse_args()


# Download taxonomy file from SILVA

# I think its better for users to input this file directly as SILVA is known to change the area that they house these files.
# We can then think about providing these files some place incase they get removed etc..
#workflow.do("curl https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_138.1.txt.gz -o [t:input/taxmap_slv_ssu_ref_138.1.txt.gz]")
#workflow.do("gunzip -fc [d:input/taxmap_slv_ssu_ref_138.1.txt.gz] > [t:input/taxmap_slv_ssu_ref_138.1.txt]")

# Add prefix of db name
dbname = Path(args.database).stem

#add the name for the generated trimmed aligned sequences
args.trimmedDatabase = os.path.join(args.output, dbname + '.pcr.align')
#add the name for the generated primer trimmed tree
args.tree = os.path.join(args.output, "region_specific.tree")
#add the name for the tree log file.
args.treelog = os.path.join(args.output, 'treelog.txt')

os.environ["OMP_NUM_THREADS"]=args.threads

## Trim database
workflow.add_task(
    "mothur -q '#set.dir(output=[args[0]]);pcr.seqs(fasta = [depends[0]], oligos=[depends[1]], pdiffs=0, rdiffs=0, keepdots=t)'",
    depends=[args.database, args.primers],
    targets=[args.trimmedDatabase],
    args=args.output,
    name="Trimming database")

## Create tree from region-specific alignment, save log file for use by pplacer

workflow.add_task(
    "FastTree -log [targets[0]] \
    -nt [depends[0]] > [targets[1]]",
    depends= [args.trimmedDatabase],
    targets=[args.treelog, args.tree],
    name="Building tree from primer-trimmed database")

## Find best thresholds

# See R script for comments
workflow.add_task(
    "src/find.cutoffs_flex_parallel.R   -d [depends[1]] -o [args[0]] -n [depends[2]] --wt1 [args[1]] --wt2 [args[2]] --threads [args[3]]",
    depends=[TrackedExecutable("src/analysis.R"), args.taxonomy, args.tree],
    targets= [args.output+"/optimal_scores.png"],
    args=[args.output, args.sweight, args.mweight, args.threads],
    name="Finding thresholds"
)


## Assign taxonomy to nodes 
workflow.add_task(
    "src/assign.node.tax_flex.R   -d [depends[1]] -o [args[0]] -n [depends[2]] --bError [args[1]] --bThreshold [args[2]]",
    depends=[TrackedExecutable("src/analysis.R"), args.taxonomy, args.tree,
             args.output+"/optimal_scores.png"],
    targets= args.output+"/resultTree_bestThresholds.RData",
    args=[args.output, args.errorRate, args.binoThreshold],
    name="Assigning taxonomy to internal nodes of ref tree"
    )

# Run the workflow
workflow.go()
