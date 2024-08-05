#!/usr/bin/python


# update parathaa code to take in one directory from step 1 with all the required files?
import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable
from pathlib import Path

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.1.1",                    #Update the version as needed
    description="Preserving Primer-Induced Taxonomic Ambiguities for Amplicons"     #Update the description as needed
    ) 

#this will take in a directory and then create all the rest of the arguments..
workflow.add_argument(
    name="treeFiles",
    desc="The output directory generated from run_tree_analysis.py",
    default=""
)

#Setting additional custom arguments for workflow - run.py
workflow.add_argument(
    name="trimmedDatabase",
    desc="",
    default=""
)
    
workflow.add_argument(
    name="trimmedTree",
    desc="Database for taxonomy [default: input/silva.seed_v138_1/silva.seed_v138_1.align]",
    default="Pre-computed-trees/Silva_v138/silva.seed_v138_1.align"
)

workflow.add_argument(
    name="treeLog",
    desc="Log file for generated tree",
    default=""
)

workflow.add_argument(
    name="query",
    desc="Reads to be taxonomically classified [default: input/test_reads_V4.fasta]",
    default="input/test_reads_V4.fasta")

workflow.add_argument(
    name="thresholds",
    desc="Rdata file that contains information on the optimal distance threshold for taxonomic level assignment",
    default=""
)

workflow.add_argument(
    name="namedTree",
    desc="Rdata file containing the formatted tree with taxonomic labels",
    default=""
)

workflow.add_argument(
    name="threads",
    desc="Number of threads to run job",
    default="1"
)

workflow.add_argument(
    name="delta",
    desc="Second thresholding value for Species that looks at neighbouring tips",
    default="0"

)

workflow.add_argument(
    name="mult",
    desc="Species threshold multiplier",
    default="0.1"
)

workflow.add_argument(
   "clean",
   desc="Clean intermediate files",
   action="store_true"
)

workflow.add_argument(
    name="deltamult",
    desc="Mutlipler for delta value if delta is not set",
    default="0.5"
)

workflow.add_argument(
    "sensitive", 
    desc="Run parathaa in sensitive mode using a non reduced optimal threshold", 
    action="store_true"
)

workflow.add_argument(
    name="genusmult",
    desc="mutlipler that can be applied to the genus level thresholding",
    default="1"
)


workflow.add_argument(
    name="minAlignLen",
    desc="The minimum length an alignment must be to pass the filter. If not set the default 20% of query length will be used",
    default="0"
)


# Parsing the workflow arguments
args = workflow.parse_args()

if(args.sensitive):
    print("Parathaa is running in sensitive mode setting delta to 0 and species threshold multipler to 1")
    args.mult=1

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
    if(len(args.treeFiles)!=0):
        db_files = [f for f in os.listdir(args.treeFiles) if f.endswith('.pcr.align')]
        for i in db_files:
            if not i.endswith('.scrap.pcr.align'):
                args.trimmedDatabase = args.treeFiles + "/" + i
        
        temp = [f for f in os.listdir(args.treeFiles) if f.endswith('.tree')]
        args.trimmedTree = args.treeFiles + "/" + temp[0]
        
        temp = [f for f in os.listdir(args.treeFiles) if f.endswith('treelog.txt')]
        args.treeLog = args.treeFiles + "/" + temp[0]

        temp = [f for f in os.listdir(args.treeFiles) if f.endswith('scores.RData')]
        args.thresholds = args.treeFiles + "/" + temp[0]

        temp = [f for f in os.listdir(args.treeFiles) if f.endswith('holds.RData')]
        args.namedTree = args.treeFiles + "/" + temp[0]


    queryName = Path(args.query).stem
    alignName = os.path.join(args.output, queryName + '.align')

    final_out = os.path.join(args.output, "taxonomic_assignments.tsv")
    alignment_report=os.path.join(args.output, queryName + '.align_report')
    ## Align query reads to trimmed seed alignment
    workflow.add_task(
        "mothur '#set.dir(output=[args[0]]);set.dir(debug=[args[0]]);align.seqs(candidate=[depends[0]], template=[depends[1]], processors=[args[1]])'",
        depends=[args.query,args.trimmedDatabase],
        targets=[alignName, alignment_report],
        args=[args.output, args.threads],
        name="Aligning queries to trimmed db"
    )

    merged = os.path.join(args.output, "merged.fasta")

    

    

    

    ## Merge Files
    workflow.add_task(
        "cat [depends[0]]  [depends[1]]  > [targets[0]] ",
        depends=[alignName, args.trimmedDatabase],
        targets=[merged],
        name="Concatenting db with queries"
    )
    
    
    ### remove poorly aligned sequences based on query length and alignment length
    poor_alignments = os.path.join(args.output, "poor_query_alignments.txt")
    
    ### grab the report which is what i actually want to do...
 
    if not(args.minAlignLen==0):
        
        workflow.add_task(
            "awk -F'\t' '{ second_col = $2 + 0 ; twelfth_col = $12 + 0 ; if (twelfth_col < 0.8 * second_col || twelfth_col > 1.2 * second_col) { print $1 } }' [depends[0]] > [targets[0]]",
            depends=alignment_report,
            targets=poor_alignments,
            name="Filtering out query sequences with alignment lengths that are not within 20% of the query sequences length"
        )
        
    else:
        workflow.add_task(
            "awk -F'\t' 'NR > 1 {twelfth_col=$12 + 0; if(twelfth_col < [args[0]]) { print $1} }' [depends[0]] > [targets[0]]",
            args=args.minAlignLen,
            depends=alignment_report,
            targets=poor_alignments,
            name="Filtering out query sequences with alignment lengths less than the minimum specificed by minAlignLen"
        )
    
    
    merged_filt = os.path.join(args.output, "merged_filt.fasta")
    ##remove poorly aligned sequences
    workflow.add_task(
        "faSomeRecords -exclude [depends[0]] [depends[1]] [targets[0]]",
        depends=[merged, poor_alignments],
        targets=merged_filt,
        name="Filtering out poor query sequences from alignment"
    )
    
    mergedSub = os.path.join(args.output, "merged_sub.fasta")
    ## Replace end gap characters ('.') with '-'
    workflow.add_task(
        "sed  '/^>/! s/\./-/g' [depends[0]] > [targets[0]]",
        depends=[merged_filt],
        targets=[mergedSub],
        name="Replacing query end gap characters"
    )
                
    args.ref = os.path.join(args.output, "ref.refpkg")

    ## Create reference package for placements
    workflow.add_task(
        "taxit create -l xxx -P [args[0]] \
        --aln-fasta [depends[0]] \
        --tree-stats [depends[1]] \
        --tree-file   [depends[2]]",
        depends=[mergedSub, args.treeLog, args.trimmedTree],
        args=[args.ref],
        targets=[args.ref],
        name="Making reference package for pplacer"
    )

    ## run pplacer
    workflow.add_task(
        "pplacer --out-dir [args[0]] -c [depends[0]] --keep-at-most 99999  \
        [depends[1]]",
        depends=[args.ref, mergedSub],
        targets=[os.path.join(args.output, "merged_sub.jplace")],
        args=[args.output],
        name="Placing queries in reference tree"
    )

    # See R script for comments
        
    ## Assign taxonomy to queries
    tax_assign_parallel = get_package_file("tax.assign_parallel", "Rscript")
    nearest_neighbour = get_package_file("nearest_neighbours_parallel", "Rscript")
    workflow.add_task(
        tax_assign_parallel+"  -j [depends[0]] -o [args[0]] -t [depends[1]] -s [depends[2]] --threads [args[1]] -d [args[2]] -m [args[3]] --util1 [args[4]] --md [args[5]] --genus_mult [args[6]]",
        depends=[os.path.join(args.output, "merged_sub.jplace"), args.namedTree, args.thresholds],
        targets=final_out,
        args=[args.output, args.threads, args.delta, args.mult, nearest_neighbour, args.deltamult, args.genusmult],
        name="Assigning taxonomy to queries"
    )

    ## Clean working directory
    if (args.clean):
        workflow.add_task(
            "rm [args[0]]/merged.fasta; rm [args[0]]/merged_sub.fasta; rm mothur.*.logfile",
            depends=final_out,
            args=[args.output],
            name="Cleaning up intermediate files"
       )

    # Run the workflow
    workflow.go()

if __name__ == '__main__':
    main()
