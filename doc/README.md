# README 

### Dependencies: 
   - R version 3.6.3
   - Python version >=3.7


##### Step 1: All the defaults values are configured in `etc/config.ini`. Please add/update the existing tasks for the private analysis workflow.

##### Usage:

```python run.py --help```

```
usage: run.py [-h] [--version] [--lines LINES] -o OUTPUT [-i INPUT]
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
## Step 2: Running the workflow 
The input and the output arguments are required for the workflow. 
```
python run.py -i input/data.tsv -o output 
```