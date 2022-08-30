import os
from glob import glob
from anadama2 import Workflow
from anadama2.tracked import TrackedExecutable

# Setting the version of the workflow and short description
workflow = Workflow(
    version="0.0.1",                    #Update the version as needed
    description="Analysis Template"     #Update the description as needed
    ) 

# Setting additional custom arguments for workflow - run.py
workflow.add_argument(
    name="lines", 
    desc="Number of lines to trim [default: 10]", 
    default="10")

workflow.add_argument(
    name="metadata", 
    desc="Metadata for performing analysis [default: input/metadata.tsv]", 
    default="input/metadata.tsv")

# Parsing the workflow arguments
args = workflow.parse_args()

#Loading the config setting
args.config = 'etc/config.ini'


# AnADAMA2 example workflow.do
workflow.do("ls /usr/bin/ | sort > [t:output/global_exe.txt]")        #Command 
workflow.do("ls $HOME/.local/bin/ | sort > [t:output/local_exe.txt]") #Command 

# Task0 sample python analysis module  - src/trim.py
workflow.add_task(
    "src/trim.py --lines [args[0]] --output [targets[0]] --input "+args.input, #Command 
    depends=[TrackedExecutable("src/trim.py")],                                #Tracking executable dependencies
    targets=args.output,                                                       #Output target directory
    args=[args.lines])                                                         #Additional arguments 


# Task1 sample python visualization module - src/plot.py
workflow.add_task(
    "src/plot.py --output [targets[0]] --input "+args.input,    #Command 
    depends=[TrackedExecutable("src/plot.py")],                 #Tracking executable dependencies
    targets=args.output)                                        #Output target directory


# Task2 sample R module  - src/analysis_example.r
workflow.add_task(
    "src/analysis.R -o [targets[0]] -d "+args.metadata,     #Command 
    depends=[TrackedExecutable("src/analysis.R")],          #Tracking executable dependencies
    targets=args.output,                                    #Output target directory
    args=[args.metadata])                                   #Additional arguments 


# Task3 add_task_group  - AnADAMA2 example to execute a task on multiple input files/dependencies
multiple_input_files = glob(os.path.join(args.output, '*.txt')) #Initializing multiple input files 
output_files = [os.path.join(args.output,'data',os.path.basename(files+"_backup")) for files in multiple_input_files]
workflow.add_task_group(
    "cp [depends[0]] [targets[0]]",                            #Command 
    depends=[multiple_input_files],   #Tracking executable dependencies
    targets=output_files)                                      #Output target directory


# private python function definition 
def remove_end_tabs_function(task):
    with open(task.targets[0].name, 'w') as file_handle_out:
        for line in open(task.depends[0].name):
            file_handle_out.write(line.rstrip() + "\n")
            
            
# Task4 add_task  - AnADAMA2 example to usage of python task function 
workflow.add_task(
    remove_end_tabs_function,                       #Calling the python function  
    depends=args.input,                             #Tracking executable dependencies
    targets=args.output+"/data/data.tsv.notabs",    #Target output
    name="remove_end_tabs")


#Task5 Add the document to the workflow
pdf_report=os.path.join(os.getcwd(),args.output,"pdfReport.pdf")
workflow.add_document(
    templates="doc/template.py",
    targets=pdf_report,
    vars={
        "introduction_text": "Demo Title"
    })

# Run the workflow
workflow.go()