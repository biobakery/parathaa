
# PARATHAA
> Preserving and Assimilating Region-specific Ambiguities in Taxonomic Hierarchical Assignments for Amplicons
 
#### Dev notes.

For input file I think for now we ignore this and then we can host the input files on another platform..
We can include a script that will download these files in a static location (whether that be zenodo or hutlab public space)

Need to make a test taxonomy assignment dataset..... (not sure the best way to do this because the test making tree script won't do the best for assignment :/ 
Hmmm... (we could just use the terrible tree and just use it as test case and leave it at that)
We could also write a seperate test script that downloads an already pre-computed tree and uses that to test assignments.


Generate pre-computed files for V1V2, V3V4, EMPV4. 

#### Introduction
> PARATHAA is a tool for taxonomic assignment while retaining information about ambiguity in placements. It is still under active development.

#### Dependencies
To run Parathaa, you should be running linux/MacOS, and you should have installed:
- Python >=3.7
- R > 4.0.0
- AnADAMA2: https://pypi.org/project/anadama2/
- mothur: https://anaconda.org/bioconda/mothur
- taxtastic: https://pypi.org/project/taxtastic/
- fasttree: http://www.microbesonline.org/fasttree/#Install
- pplacer: https://anaconda.org/bioconda/pplacer


####  Cloning Repo
  

- Clone your recently created repository in your local development environment either using:
```
git clone https://github.com/biobakery/parathaa
```
or using the "**Clone or Download**" button.


##### Demo Run:
```
python run.py --primers input/V4V5.oligos --database input/silva.seed_v138_1/silva.seed_v138_1.align --query input/ASVs.fasta --output output 
```
```python run.py --help```

```
usage: run.py [-h] [--version] [--primers PRIMERS]
[--database DATABASE] [--query QUERY] -o OUTPUT [-i INPUT]
[--config CONFIG] [--local-jobs JOBS] [--grid-jobs GRID_JOBS]
[--grid GRID] [--grid-partition GRID_PARTITION]
[--grid-benchmark {on,off}] [--grid-options GRID_OPTIONS]
[--grid-environment GRID_ENVIRONMENT]
[--grid-scratch GRID_SCRATCH] [--dry-run] [--skip-nothing]
[--quit-early] [--until-task UNTIL_TASK]
[--exclude-task EXCLUDE_TASK] [--target TARGET]
[--exclude-target EXCLUDE_TARGET]
[--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
```

