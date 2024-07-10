#!/usr/bin/python

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
    desc="File with primer oligos [default: none]",
    default="none")
    
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

workflow.add_argument(
    name="clean",
    desc="Clean out intermediate files [default: true]",
    default="true"
)

# Parsing the workflow arguments
args = workflow.parse_args()
def get_package_file(basename, type="template"):
    """ Get the full path to a file included in the installed python package.

        Args:
            basename (string) : The basename of the file
            type (string) : The type of file to find (template or Rscript)

        Returns: 
            string : The full path to the file

    """

    if type == "template":
        subfolder = "document_templates"
        extension = ".pmd"
    elif type == "image":
        subfolder = os.path.join(os.path.pardir,"images")
    else:
        subfolder = "utility"
        extension = ".R"

    # get all of the templates in this folder
    package_install_folder=os.path.join(os.path.dirname(os.path.realpath(__file__)), subfolder)

    # return the template with the name
    if type != "image":
        found_files=list(filter(lambda file: file.endswith(extension),os.listdir(package_install_folder)))
        matching_file=list(filter(lambda file: file.startswith(basename+extension), found_files))
    else:
        found_files=list(filter(lambda file: os.path.isfile(os.path.join(package_install_folder,file)),os.listdir(package_install_folder)))
        matching_file=list(filter(lambda file: basename.lower() in file.lower(), found_files))

    if matching_file:
        matching_file=os.path.join(package_install_folder,matching_file[0])
    else:
        matching_file=""

    return matching_file


def main():
    # Add prefix of db name
    dbname = Path(args.database).stem

    #add the name for the generated trimmed aligned sequences
    args.trimmedDatabase = os.path.join(args.output, dbname + '.pcr.align')
    #add the name for the generated primer trimmed tree
    args.tree = os.path.join(args.output, "region_specific.tree")
    #add the name for the tree log file.
    args.treelog = os.path.join(args.output, 'treelog.txt')

    final_out = os.path.join(args.output, "resultTree_bestThresholds.RData")

    #set an environmental variable used to determine the number of threads
    os.environ["OMP_NUM_THREADS"]=args.threads

    ## Trim the aligned 16S database to the primer inputs
    if(args.primers!="none"):
        workflow.add_task(
            "mothur -q '#set.dir(output=[args[0]]);pcr.seqs(fasta = [depends[0]], oligos=[depends[1]], pdiffs=3, rdiffs=3, keepdots=t, checkorient=false, processors=[args[1]])'",
            depends=[args.database, args.primers],
            targets=[args.trimmedDatabase],
            args=[args.output, args.threads],
            name="Trimming database")
    else:
        workflow.add_task(
            "cp [args[0]] [targets[0]]",
            args=args.database,
            targets=args.trimmedDatabase
        )

    ## Create tree from region-specific alignment, save log file for use by pplacer
    workflow.add_task(
        "FastTree -log [targets[0]] \
        -nt -gtr [depends[0]] > [targets[1]]",
        depends= [args.trimmedDatabase],
        targets=[args.treelog, args.tree],
        name="Building tree from primer-trimmed database")


    # Identify the optimal distance thresholdings for taxonomic assignment
    find_cutoffs_flex_parallel = get_package_file("find.cutoffs_flex_parallel", "Rscript")
    spec_dev = get_package_file("SILVA.species.editor.dev", "Rscript")
    calc_error = get_package_file("calc.error.scores", "Rscript")
    workflow.add_task(
        find_cutoffs_flex_parallel+" -d [depends[0]] -o [args[0]] -n [depends[1]] --wt1 [args[1]] --wt2 [args[2]] --threads [args[3]] --util1 [args[4]] --util2 [args[5]]",
        depends=[args.taxonomy, args.tree],
        targets= [args.output+"/optimal_scores.png"],
        args=[args.output, args.sweight, args.mweight, args.threads, spec_dev, calc_error],
        name="Finding thresholds"
    )


    ## Assign taxonomy to the internal nodes of the phylogenetic tree
    assign_node_tax_flex = get_package_file("assign.node.tax_flex", "Rscript")
    single_tax = get_package_file("single.tax", "Rscript")
    workflow.add_task(
        assign_node_tax_flex+" -d [depends[0]] -o [args[0]] -n [depends[1]] --bError [args[1]] --bThreshold [args[2]] --util1 [args[3]] --util2 [args[4]]",
        depends=[args.taxonomy, args.tree,
                args.output+"/optimal_scores.png"],
        targets=final_out,
        args=[args.output, args.errorRate, args.binoThreshold, spec_dev, single_tax],
        name="Assigning taxonomy to internal nodes of ref tree"
        )

    if args.clean == "true":
        workflow.add_task(
            "rm [args[0]]/internal_node_stats.RData; rm [args[0]]/*.bad.accnos; rm [args[0]]/*.pcr.8mer; rm [args[0]]/*.scrap.pcr.align; rm mothur.*.logfile",
            depends=final_out,
            args=args.output,
            name="Cleaning intermediate files"
        )

    # Run the workflow
    workflow.go()

if __name__ == '__main__':
    main()
