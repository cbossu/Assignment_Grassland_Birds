#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=pcangsd
#################
#a file for job output, you can check job progress
#SBATCH --output=pcangsd.%j.out
#################
# a file for errors from the job
#SBATCH --error=pcangsd.%j.err
#################
#time you think you need; default is one hour
#in minutes in this case
#SBATCH -t 8:00:00
#################
#partition and account
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --account=csu88_alpine2
#################
#number of threads you want to use and memory
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
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
conda activate bcftools
plink --bfile input/feems/LOSH.top_pairwiseFST.second.exclLD --aec --keep input/Breeding_ind.plink.txt \
  --out input/feems/LOSH.top_pairwiseFST.second.exclLD.Breed_only --make-bed

conda activate pcangsd
pcangsd -p input/feems/LOSH.top_pairwiseFST.second.exclLD.Breed_only -o output/feems/LOSH.top_pairwiseFST.second.exclLD.Breed_only.k8.bed.pcangsd -t 4

