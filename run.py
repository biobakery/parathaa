import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.1.1",                    #Update the version as needed
    description="Preserving Primer-Induced Taxonomic Ambiguities for Amplicons"     #Update the description as needed
    ) 

# Setting additional custom arguments for workflow - run.py
workflow.add_argument(
    name="primers",
    desc="File with primer oligos [default: input/EMPV4.oligos]",
    default="input/EMPV4.oligos")
    
workflow.add_argument(
    name="database",
    desc="Database for taxonomy [default: input/silva.seed_v138_1/silva.seed_v138_1.align]",
    default="input/silva.seed_v138_1/silva.seed_v138_1.align")

workflow.add_argument(
    name="query",
    desc="Reads to be taxonomically classified [default: input/test_reads_V4.fasta]",
    default="input/test_reads_V4.fasta")

# Parsing the workflow arguments
args = workflow.parse_args()

#Loading the config setting
args.config = 'etc/config.ini'


# Download taxonomy file from SILVA

# *Hmmm this is probably better to be in the config file encase this link changes?
# Also I wonder if you want to include downloading this as an optional step incase people want to use different taxonomy files?
workflow.do("curl https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/taxmap_slv_ssu_ref_138.1.txt.gz -o [t:input/taxmap_slv_ssu_ref_138.1.txt.gz]")
workflow.do("gunzip -fk [d:input/taxmap_slv_ssu_ref_138.1.txt.gz] > [t:input/taxmap_slv_ssu_ref_138.1.txt]")

# Adding name of trimmed Seqs
dbname = args.database

# interesting that the args.trimmedDatabase is placed in the input folder for that database... is there a way we can fix this?
# perhaps running the mothur set.dir function before running the #pcr seqs command....

args.trimmedDatabase = dbname[:dbname.find(".align")] + '.pcr.align'
args.tree = os.path.join(args.output, "region_specific.tree")
args.treelog = os.path.join(args.output, 'treelog.txt')

## Trim database

# I wonder for this step if there could be a log file showing how many sequences in the database passed this step
# from what I understand the database has to be a pre-aligned for input into parathaa?

workflow.add_task(
    "mothur '#pcr.seqs(fasta = [depends[0]], oligos=[depends[1]], pdiffs=0, rdiffs=0, keepdots=f)'",
    depends=[args.database, args.primers],
    targets=args.trimmedDatabase,
    name="Trimming database")

## Create tree from region-specific alignment, save log file for use by pplacer

# how does this step deal with duplicate reads that are generated after primer trimming?

workflow.add_task(
    "fasttree -gtr -log [targets[0]] \
    -nt [depends[0]] > [targets[1]]",
    depends= args.trimmedDatabase,
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

queryName = args.query
alignName = queryName[:queryName.find(".fasta")] + '.align'



## Align query reads to trimmed seed alignment
workflow.add_task(
    "mothur '#align.seqs(candidate=[depends[0]], template=[depends[1]])'",
    depends=[args.query,args.trimmedDatabase],
    targets=alignName,
    name="Aligning queries to trimmed db"
    )

merged = os.path.join(args.output, "merged.fasta")
mergedSub = merged[:merged.find(".fasta")] + '.sub.fasta'

## Merge Files
workflow.add_task(
    "cat [depends[0]]  [depends[1]]  > [targets[0]] ",
    depends=[alignName, args.trimmedDatabase],
    targets=[merged],
    name="Concatenting db with queries"
    )

## Replace end gap characters ('.') with '-'
workflow.add_task(
    "sed  '/^>/! s/\./-/g' [depends[0]] > [targets[0]]",
    depends=[merged],
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
    depends=[mergedSub, args.treelog, args.tree],
    args=[args.ref],
    targets=[args.ref],
    name="Making reference package for pplacer"
    )

## run pplacer
workflow.add_task(
    "pplacer --out-dir [args[0]] -c [depends[0]]  \
    [depends[1]]",
    depends=[args.ref, mergedSub],
    targets=[os.path.join(args.output, "merged.sub.jplace")],
    args=[args.output],
    name="Placing queries in reference tree"
    )

# See R script for comments
    
## Assign taxonomy to queries
workflow.add_task(
    "src/tax.assign.R  -j [depends[0]] -o [args[0]] -t [depends[1]]",
    depends=[os.path.join(args.output, "merged.sub.jplace"), args.output+"/resultTree_bestThresholds.RData"],
    targets=[os.path.join(args.output, "taxonomic_assignments.tsv")],
    args=[args.output],
    name="Assigning taxonomy to queries"
    )

# Run the workflow
workflow.go()
