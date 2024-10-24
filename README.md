# Chienlab-breseq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**chienlab-breseq** is a simple nextflow wrapper of the [Breseq](https://github.com/barricklab/breseq) pipeline meant for parallel execution on high performance clusters. 


## Installation

You will need to install [`Nextflow`](https://www.nextflow.io/) (version 21.10.3+).

```
Usage:
nextflow run baldikacti/chienlab-breseq --project [dir] --gbk [file] -profile conda [other_options]

Mandatory arguments:
  --project [directory]           Path to directory containing FastQ files.
  --gbk [file]                    Path to GBK file containing reference genome.
  -profile [str]                  Configuration profile to use.
                                  Available: conda

Other options:
  --gbk2 [file]                   Path to a second gbk formatted reference file.
  --max_memory ['32.GB']          Maximum available memory in the system.
  --max_cpus [int]                Maximum available cpu's in the system.
  --max_time ['10.h']             Maximum time requested time for the pipeline.
  -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
  -log [file]                     Path to optional log file
```

## Example Usage

1. Put your `.fastq` files in a directory. In this example we will put out files in `data/test` directory.

NOTE: DO NOT change the name of your fastq files.

**Example folder structure:**

Here we have two samples `test` and `test2` in the `data/` folder.

```
data/
├── test
│   ├── test_R1.fastq.gz
│   └── test_R2.fastq.gz
└── test2
    ├── test2_R1.fastq.gz
    └── test2_R2.fastq.gz
```

2. Open terminal in where you want to run the nextflow pipeline.

3. Setup a `bash` script as below. 

With single reference file.
**main.sh**
```bash
#!/usr/bin/bash

set -e 

nextflow run baldikacti/chienlab-breseq \
    --project data/test1  \
    --gbk references/genome1.gbk \
    -profile conda
```

With two reference files.
**main.sh**
```bash
#!/usr/bin/bash

set -e 

nextflow run baldikacti/chienlab-breseq \
    --project data/test1  \
    --gbk references/genome1.gbk \
    --gbk2 references/genome2.gbk \
    -profile conda
```

4. Run

```bash
bash main.sh
```

## Slurm Example

This example uses settings for the [Unity HPC](https://unity.rc.umass.edu/) using the **slurm** scheduler. Adjust settings for different clusters and schedulers.

**main_slurm.sh**
```bash
#!/usr/bin/bash
#SBATCH --job-name=chienlab-breseq             # Job name
#SBATCH --partition=cpu,cpu-preempt            # Partition (queue) name
#SBATCH -c 24                                  # Number of CPUs
#SBATCH --nodes=1                              # Number of nodes
#SBATCH --mem-per-cpu=4gb                      # Job memory request
#SBATCH --time=06:00:00                        # Time limit hrs:min:sec
#SBATCH --output=logs/chienlab-breseq_%j.log   # Standard output and error log

set -e

module load nextflow/24.04.3 conda/latest

nextflow run baldikacti/chienlab-breseq \
    --project data/test1  \
    --gbk references/NC_011916.gbk \
    -profile conda
```

And submit the job with
```bash
sbatch main.sh
```