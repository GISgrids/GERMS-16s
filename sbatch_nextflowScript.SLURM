#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=li_zhihui@gis.a-star.edu.sg
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --error=ferrero_16s.%J.err
#SBATCH --output=ferrero_16s.%J.out

# Define location of your own location where the 16s location is
S3_DIR=$HOME/GERMS_16s_resource

# Download data to scratch first and change directory
cd /scratch/users/astar/gis/lizh/Projects/GERMS-16s

# TEST data
nextflow run ${S3_DIR}/GERMS_16s_ORIGINAL_ACRC_v1.nf -resume -profile acrc --reads "test/*/*_R{1,2}_unaligned_*.fastq.gz" --outdir 'test_16s_HOME' --s3_dir ${S3_DIR}

# nextflow run $HOME/GERMS_16s_ORIGINAL_ACRC_v1.nf -resume -profile acrc --reads "fastq_ferrero/*/*_R{1,2}_unaligned_*.fastq.gz" --outdir 'results_16s_Ferroro'
