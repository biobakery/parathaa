
# PARATHAA
> Preserving and Assimilating Region-specific Ambiguities in Taxonomic Hierarchical Assignments for Amplicons
 
#### Introduction
> PARATHAA is a tool for taxonomic assignment while retaining information about ambiguity in placements. It is still under active development.

#### Dependencies
To run Parathaa, you should be running linux/MacOS, and you should have installed:
- Python >=3.7
- R > 4.0.0
- AnADAMA2: https://pypi.org/project/anadama2/
- taxtastic: https://pypi.org/project/taxtastic/
- pplacer: https://anaconda.org/bioconda/pplacer
- fasttree: http://www.microbesonline.org/fasttree/#Install


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

