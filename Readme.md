
# PARATHAA
> Preserving and Assimilating Region-specific Ambiguities in Taxonomic Hierarchical Assignments for Amplicons
 
#### Introduction
> PARATHAA is a tool for taxonomic assignment while retaining information about ambiguity in placements. It is still under active development.

  
####  Cloning Repo
  

- Clone your recently created repository in your local development environment either using:
```
git clone https://github.com/biobakery/parathaa
```
or using the "**Clone or Download**" button.


##### Demo Run:
```
python run.py --input input/data.tsv --output output --lines 10 --metadata input/metadata.tsv 
```
```python run.py --help```

```
usage: run.py [-h] [--version] [--lines LINES]
[--sample-metadata SAMPLE_METADATA] -o OUTPUT [-i INPUT]
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

