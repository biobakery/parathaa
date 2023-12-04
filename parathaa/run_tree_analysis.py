#!/usr/bin/python3

import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable
from pathlib import Path


# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.2.1",
    description=("Preserving Primer-Induced Taxonomic Ambiguities for Amplicons: This script given an aligned 16S file, an accompanying taxonomy file "
                 " and a primer pair will output a phylogenetic tree with taxonomically labeled internal nodes based on the amplified region using optimal "
                 "distance thresholds for taxonomic assignment." 
    ) 
)

# Setting additional custom arguments for workflow - run_tree_analysis
# Default parameters will run the test case shown in the tutorial.

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
    desc="File that contains taxonomy for the input database [defualt: input/taxmap_slv_ssu_ref_138.1.txt]",
    default="input/taxmap_slv_ssu_ref_138.1.txt"
)

workflow.add_argument(
    name="sweight",
    desc="penalty weight for over-splitting errors when identifying optimal distance thresholds [default: 1]",
    default=1
)

workflow.add_argument(
    name="mweight",
    desc="penalty weight for over-merging errors when identifying optimal distance thresholds [default: 1]",
    default=1
)

workflow.add_argument(
    name="errorRate",
    desc="The assumed error rate for the binomial model for ambigious taxonomic assignment [default: 0.05]",
    default=0.05
)

workflow.add_argument(
    name="binoThreshold",
    desc="Critial p-value (threshold) for the binomial error model for ambigious taxonomic assignment [default: 0.20]",
    default=0.20
)


workflow.add_argument(
    name="threads",
    desc="Number of threads to run multi-threaded processes [default: 1]",
    default="1"

)

# This needs a better descripition from Meg
workflow.add_argument(
    name="namePar",
    desc="Testing parameter for Silva name editor [default: spNames4]",
    default="spNames4"
)

# Parsing the workflow arguments
args = workflow.parse_args()

def main():
    # Add prefix of db name
    dbname = Path(args.database).stem

    #add the name for the generated trimmed aligned sequences
    args.trimmedDatabase = os.path.join(args.output, dbname + '.pcr.align')
    #add the name for the generated primer trimmed tree
    args.tree = os.path.join(args.output, "region_specific.tree")
    #add the name for the tree log file.
    args.treelog = os.path.join(args.output, 'treelog.txt')

    #set an environmental variable used to determine the number of threads
    os.environ["OMP_NUM_THREADS"]=args.threads

    ## Trim the aligned 16S database to the primer inputs
    workflow.add_task(
        "mothur -q '#set.dir(output=[args[0]]);pcr.seqs(fasta = [depends[0]], oligos=[depends[1]], pdiffs=3, rdiffs=3, keepdots=t)'",
        depends=[args.database, args.primers],
        targets=[args.trimmedDatabase],
        args=args.output,
        name="Trimming database")

    ## Create tree from region-specific alignment, save log file for use by pplacer
    workflow.add_task(
        "FastTree -log [targets[0]] \
        -nt -gtr [depends[0]] > [targets[1]]",
        depends= [args.trimmedDatabase],
        targets=[args.treelog, args.tree],
        name="Building tree from primer-trimmed database")


    # Identify the optimal distance thresholdings for taxonomic assignment
    workflow.add_task(
        "utility/find.cutoffs_flex_parallel.R   -d [depends[0]] -o [args[0]] -n [depends[1]] --wt1 [args[1]] --wt2 [args[2]] --threads [args[3]] --name [args[4]]",
        depends=[args.taxonomy, args.tree],
        targets= [args.output+"/optimal_scores.png"],
        args=[args.output, args.sweight, args.mweight, args.threads, args.namePar],
        name="Finding thresholds"
    )


    ## Assign taxonomy to the internal nodes of the phylogenetic tree
    workflow.add_task(
        "utility/assign.node.tax_flex.R   -d [depends[0]] -o [args[0]] -n [depends[1]] --bError [args[1]] --bThreshold [args[2]]",
        depends=[args.taxonomy, args.tree,
                args.output+"/optimal_scores.png"],
        targets= args.output+"/resultTree_bestThresholds.RData",
        args=[args.output, args.errorRate, args.binoThreshold],
        name="Assigning taxonomy to internal nodes of ref tree"
        )


    # Run the workflow
    workflow.go()

if __name__ == '__main__':
    main()