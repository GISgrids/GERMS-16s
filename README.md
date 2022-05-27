# GERMS-16s
This is a repo for the original 16S pipeline developed with in the Genome Institute of Singapore (http://www.a-star.edu.sg/gis), designed for Illumina shotgun sequencing of 16S rRNA amplicon sequences ([Ong et al., 2013, PMID 23579286](http://www.ncbi.nlm.nih.gov/pubmed/23579286)).

This current repo is a remake of pipeliune with Nextflow in DSL1 and is intended for running in SLURM 

The original Snakemake version of the pipeline can be found [CSB5/GERMS_16S_pipeline] (https://github.com/CSB5/GERMS_16S_pipeline)

## Setup 
Contact author for S3 bucket which contains the database and script in the correct location.
OR
DIY
1) install nextflow 
2) install conda - preferably miniconda
3) Prepare Greengene 
    a) Create and go into folder database/greengenes/
        `mkdir -p database/greengenes/`
        `cd database/greengenes/`
    b) Download [gg_13_5_tax_with_hOTUs_99_reps.fasta](https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_tax_with_hOTUs_99_reps.fasta)
    c) Create index for [Graphmap2](https://github.com/lbcb-sci/graphmap2)
4) Prepare [emirge](https://github.com/csmiller/EMIRGE) 
    a) Create and go into folder database/ssu/
       `mkdir -p database/greengenes/`
       `cd database/greengenes/`
    b) `emirge_makedb.py -p 16 -r 138.1 -m 800 -M 3500`
  
## Running pipeline in SLURM cluster
Edit the file 'sbatch_nextflowScript.SLURM' to change directory (cd) into the work directory where this repo has been cloned to
then 
`sbatch sbatch_nextflowScript.SLURM`
