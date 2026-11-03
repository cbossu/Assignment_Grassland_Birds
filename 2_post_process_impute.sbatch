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
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --account=csu88_alpine2
#SBATCH --ntasks-per-node=4
#SBATCH --mem=11G
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
conda activate feems_e


python scripts/filter_gprobs_by_r2.py \
        --gprobs output/imputation/losh.impute3.combined.subsetted.beagle.gz.gprobs.gz \
        --r2 output/imputation/losh.impute3.combined.subsetted.beagle.gz.r2 \
        --min-r2 0.6 \
        --out-gprobs output/imputation/losh.impute3.filtered.gprobs.gz \
        --hard-calls output/imputation/losh.impute3.filtered.geno.tsv.gz \
        --vcf output/imputation/losh.impute3.filtered.vcf \
        --min-call-prob 0.8 \
        --sample-names input/samples.txt          # optional, one ID per column, in order
        
conda activate bcftools
mkdir input/pairwiseFST/

## Remove one duplicated individual
bcftools view -Oz -o input/pairwiseFST/losh.impute3.filtered.vcf.gz  -s ^21N02103b --threads 4 output/imputation/losh.impute3.filtered.vcf
tabix -p vcf input/pairwiseFST/losh.impute3.filtered.vcf.gz
rm output/imputation/losh.impute3.filtered.vcf

## Split population into train (90%) and test (10%) versions
conda activate R-GEA2
Rscript 3_split_train_test_fst.R

## Pairwise fst between groups in pop folder
conda activate bcftools
bash output/fst_splits/rep01/run_fst.sh

# Select top outlier variatns amongst all pairwise-FST analyses
conda activate R-GEA2
Rscript 4_select_top_fst.R output/fst_splits/rep01


## Subset to outlier loci
conda activate bcftools
 bcftools view -T output/fst_splits/rep01/top_fst_positions.train.txt input/pairwiseFST/losh.impute3.filtered.vcf.gz -Oz -o output/fst_splits/rep01/LOSH_topfst.vcf.gz

## Convert to plink for feems, and remove linked loci
bash 5_vcf_to_plink.sh output/fst_splits/rep01/LOSH_topfst.vcf.gz input/feems/LOSH.feems_meta.csv input/feems/LOSH.top_pairwiseFST.second
