#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=impute_vcf
#################
#a file for job output, you can check job progress
#SBATCH --output=impute.%j.out
#################
# a file for errors from the job
#SBATCH --error=impute.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 24:00:00
#################
#number of threads you want to use
#SBATCH --cpus-per-task 4
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=cbossu@rams.colostate.edu
#################
#now run normal batch commands
##################
#echo commands to stdout

set -x

source ~/.bashrc
conda activate imputation


python scripts/filter_gprobs_by_r2.py \
        --gprobs output/imputation/losh.impute3.combined.subsetted.beagle.gz.gprobs.gz \
        --r2 output/imputation/losh.impute3.combined.subsetted.beagle.gz.r2 \
        --min-r2 0.6 \
        --out-gprobs output/imputation/losh.impute3.filtered.gprobs.gz \
        --hard-calls output/imputation/losh.impute3.filtered.geno.tsv.gz \
        --vcf output/imputation/losh.impute3.filtered.vcf.gz \
        --min-call-prob 0.8 \
        --sample-names input/samples.txt          # optional, one ID per column, in order
        
## Annotate the vcf file
conda activate bcftools
bcftools annotate --set-id '%CHROM\:%POS'  output/imputation/losh.impute3.filtered.vcf.gz \
  -Oz -o input/pairwiseFST/losh.impute3.filtered.annotated.vcf.gz
