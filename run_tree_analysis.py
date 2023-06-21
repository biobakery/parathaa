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

# Parsing the workflow arguments
args = workflow.parse_args()

#Loading the config setting
args.config = 'etc/config.ini'


# Download taxonomy file from SILVA

# I think its better for users to input this file directly as SILVA is known to change the area that they house these files.
# We can then think about providing these files some place incase they get removed etc..
#workflow.do("curl https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_138.1.txt.gz -o [t:input/taxmap_slv_ssu_ref_138.1.txt.gz]")
#workflow.do("gunzip -fc [d:input/taxmap_slv_ssu_ref_138.1.txt.gz] > [t:input/taxmap_slv_ssu_ref_138.1.txt]")

# Add prefix of db name
dbname = Path(args.database).stem
print(dbname)

#add the name for the generated trimmed aligned sequences
args.trimmedDatabase = os.path.join(args.output, dbname + '.pcr.align')
#add the name for the generated primer trimmed tree
args.tree = os.path.join(args.output, "region_specific.tree")
#add the name for the tree log file.
args.treelog = os.path.join(args.output, 'treelog.txt')


## Trim database
workflow.add_task(
    "mothur '#set.dir(output=[args[0]]);pcr.seqs(fasta = [depends[0]], oligos=[depends[1]], pdiffs=0, rdiffs=0, keepdots=f)'",
    depends=[args.database, args.primers],
    targets=[args.trimmedDatabase],
    args=args.output,
    name="Trimming database")

## Create tree from region-specific alignment, save log file for use by pplacer

# how does this step deal with duplicate reads that are generated after primer trimming?

workflow.add_task(
    "FastTree -gtr -log [targets[0]] \
    -nt [depends[0]] > [targets[1]]",
    depends= [args.trimmedDatabase],
    targets=[args.treelog, args.tree],
    name="Building tree from primer-trimmed database")

## Find best thresholds

# See R script for comments
# Does this try and find optimal cutoffs for negative binomial model?
workflow.add_task(
    "src/find.cutoffs.R   -d [depends[1]] -o [args[0]] -n [depends[2]]",
    depends=[TrackedExecutable("src/analysis.R"), "input/taxmap_slv_ssu_ref_138.1.txt", args.tree],
    targets= [args.output+"/optimal_scores.png"],
    args=args.output,
    name="Finding thresholds"
)


# See R script for comments

# did my best to go through this script but still a bit unclear on some parts of it!

## Assign taxonomy to nodes 
## it was a bit unclear to me but is the point of this to determine the taxonomy of internal nodes 
## based on the new tree? so we can use that to assign taxonomy to the newly placed tips?
workflow.add_task(
    "src/assign.node.tax.R   -d [depends[1]] -o [args[0]] -n [depends[2]]",
    depends=[TrackedExecutable("src/analysis.R"), "input/taxmap_slv_ssu_ref_138.1.txt", args.tree,
             args.output+"/optimal_scores.png"],
    targets= args.output+"/resultTree_bestThresholds.RData",
    args=args.output,
    name="Assigning taxonomy to internal nodes of ref tree"
    )

# Run the workflow
workflow.go()
