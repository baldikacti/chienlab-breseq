# Breseq Protocol

1. Upload your data to `data/` folder.

Your fastq files needs to be directly in the folder you upload to `data/`.

DO NOT change the name of your fastq files.


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

2. Open terminal in the breseq folder. 

3. Run the following command in the terminal
Replace `data/test` with the path to your data folder.

```
sbatch main_slurm.sh data/test
```