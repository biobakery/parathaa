import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.0.1",                    #Update the version as needed
    description="Preserving Primer-Induced Taxonomic Ambiguities for Amplicons"     #Update the description as needed
    ) 

# Setting additional custom arguments for workflow - run.py
workflow.add_argument(
    name="primers",
    desc="File with primer oligos [default: input/EMPV4.oligos]",
    default="input/EMPV4.oligos")

workflow.add_argument(
    name="inputSeqs",
    desc="Optional sequences to be assigned taxonomy",
    default=0)
    
workflow.add_argument(
    name="database",
    desc="Database for taxonomy [default: input/silva.seed_v138_1/silva.seed_v138_1.align]",
    default="input/silva.seed_v138_1/silva.seed_v138_1.align")

# Parsing the workflow arguments
args = workflow.parse_args()

#Loading the config setting
args.config = 'etc/config.ini'


# AnADAMA2 example workflow.do
#workflow.do("ls /usr/bin/ | sort > [t:output/global_exe.txt]")        #Command 
#workflow.do("ls $HOME/.local/bin/ | sort > [t:output/local_exe.txt]") #Command 

# Task0 sample python analysis module  - src/trim.py
#workflow.add_task(
#    "src/trim.py --lines [args[0]] --output [targets[0]] --input "+args.input, #Command 
#    depends=[TrackedExecutable("src/trim.py")],                                #Tracking executable dependencies
#    targets=args.output,                                                       #Output target directory
#    args=[args.lines])                                                         #Additional arguments 

# Adding name of trimmed Seqs
dbname = args.database
args.trimmedDatabase = dbname[:dbname.find(".align")] + '.pcr.align'
args.tree = "region_specific.tree"

## Trim database
workflow.add_task("/Applications/mothur/mothur '#pcr.seqs(fasta = [depends[0]], oligos=[depends[1]], pdiffs=0, rdiffs=0)'",
                  depends=[args.database, args.primers],
                  targets=args.trimmedDatabase,
                  name="trim database")

## Create tree from region-specific alignment, save log file for use by pplacer
workflow.add_task("/Users/mis696/anaconda3/lib/python3.8/site-packages/apples/tools/FastTree-darwin -gtr -log SILVA_seed_V4_log.txt \
	   -nt [depends[0]] > [targets[0]]",
                  depends=args.trimmedDatabase,
                  targets=args.tree)


# Task2 sample R module  - src/analysis_example.r
#workflow.add_task(
#    "src/analysis.R -o [targets[0]] -d "+args.metadata,     #Command 
#    depends=[TrackedExecutable("src/analysis.R")],          #Tracking executable dependencies
#    targets=args.output,                                    #Output target directory
#    args=[args.metadata])                                   #Additional arguments 


# Task3 add_task_group  - AnADAMA2 example to execute a task on multiple input files/dependencies
#multiple_input_files = glob(os.path.join(args.output, '*.txt')) #Initializing multiple input files 
#output_files = [os.path.join(args.output,'data',os.path.basename(files+"_backup")) for files in multiple_input_files]
#workflow.add_task_group(
#    "cp [depends[0]] [targets[0]]",                            #Command 
#    depends=[multiple_input_files],   #Tracking executable dependencies
#    targets=output_files)                                      #Output target directory


# private python function definition 
#def remove_end_tabs_function(task):
#    with open(task.targets[0].name, 'w') as file_handle_out:
#        for line in open(task.depends[0].name):
#            file_handle_out.write(line.rstrip() + "\n")
            
            
# Task4 add_task  - AnADAMA2 example to usage of python task function 
#workflow.add_task(
#    remove_end_tabs_function,                       #Calling the python function  
#    depends=args.input,                             #Tracking executable dependencies
#    targets=args.output+"/data/data.tsv.notabs",    #Target output
#    name="remove_end_tabs")


# Run the workflow
workflow.go()
