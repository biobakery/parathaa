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

# Setting additional custom arguments for workflow - run.py
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

# Parsing the workflow arguments
args = workflow.parse_args()

queryName = Path(args.query).stem
alignName = os.path.join(args.output, queryName + '.align')

## Align query reads to trimmed seed alignment
workflow.add_task(
    "mothur '#set.dir(output=[args[0]]);set.dir(debug=[args[0]]);align.seqs(candidate=[depends[0]], template=[depends[1]])'",
    depends=[args.query,args.trimmedDatabase],
    targets=alignName,
    args=args.output,
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
mergedSub = os.path.join(args.output, "merged_sub.fasta")
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
    depends=[mergedSub, args.treeLog, args.trimmedTree],
    args=[args.ref],
    targets=[args.ref],
    name="Making reference package for pplacer"
)

## run pplacer
workflow.add_task(
    "pplacer --out-dir [args[0]] -c [depends[0]]  \
    [depends[1]]",
    depends=[args.ref, mergedSub],
    targets=[os.path.join(args.output, "merged_sub.jplace")],
    args=[args.output],
    name="Placing queries in reference tree"
)

# See R script for comments
    
## Assign taxonomy to queries
workflow.add_task(
    "src/tax.assign.R  -j [depends[0]] -o [args[0]] -t [depends[1]] -s [depends[2]]",
    depends=[os.path.join(args.output, "merged_sub.jplace"), args.namedTree, args.thresholds],
    targets=[os.path.join(args.output, "taxonomic_assignments.tsv")],
    args=[args.output],
    name="Assigning taxonomy to queries"
)

# Run the workflow
workflow.go()
