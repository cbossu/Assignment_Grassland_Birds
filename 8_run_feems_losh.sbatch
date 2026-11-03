#!/bin/bash
#SBATCH --job-name=feems_losh
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --account=csu88_alpine2
#SBATCH --ntasks-per-node=8
#SBATCH --mem=29G
#SBATCH --time=24:00:00
#SBATCH --output=feems_losh_%j.log
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

# CV (leave-one-deme-out by default) and --loo each refit FEEMS once per sampled deme,
# so runtime scales with the number of sampled demes. Wolves (72 demes, ~17k SNPs):
# ~9 min for 5-fold CV over the 10x5 lambda grid on 2 cores.

source ~/.bashrc
conda activate feems_e
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
export OPENBLAS_NUM_THREADS=$SLURM_NTASKS_PER_NODE

## Can run on acompile on alpine. If so, use the code below, if run on alpine, comment out the export functions below
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4

GRID=/projects/cbossu@colostate.edu/feems/feems/data/grid_250.shp   # ~220 km spacing; grid_100 if breeding sites are dense
GRID=/projects/cbossu@colostate.edu/feems/feems/data/grid_100.shp   # grid_100 if breeding sites are dense (test- LOSH has a lot of breeding samples)

python scripts/feems_assign_nonbreeding.py \
  --bfile input/feems/LOSH.top_pairwiseFST.second.exclLD \
  --meta  input/feems/LOSH.feems_meta.csv \
  --grid  $GRID \
  --outer input/feems/LOSH_outer.txt \
  --buffer 2 \
  --cv_folds 10 \
  --method point_mu \
  --loo  --loo_ind  --ind_folds 10 \
  --outdir output/feems/feems_LOSH_fst_outer

# Second run without the --loo or --loo_ind, using the T
python scripts/feems_assign_nonbreeding.py \
  --bfile input/feems/LOSH.top_pairwiseFST.second.exclLD \
  --meta  input/feems/LOSH.feems_meta.csv \
  --grid  $GRID \
  --outer input/feems/LOSH_outer.txt \
  --buffer 2 \
  --cv_folds 10 \
  --method point_mu \
  --temp 251 \
  --outdir output/feems/feems_LOSH_fst_outer_T251


conda activate R-GEA2
Rscript scripts/plot_feems.R output/feems/feems_LOSH_fst_outer_T251/ output/feems/feems_LOSH_fst_outer/