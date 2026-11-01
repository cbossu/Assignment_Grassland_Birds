#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=impute
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
conda init
conda activate imputation

PROJ_DIR="/home/cbossu/programs/beagle_v3.3/"
java -Djava.io.tmpdir="/home/cbossu/LOSH/OriGen/tmp_dir" -Xmx14000m -jar "$PROJ_DIR"/beagle.jar \
   like=input/combined.subsetted.beagle.gz out=output/imputation/losh.impute3
   
