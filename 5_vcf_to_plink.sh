#!/usr/bin/env bash
# Beagle-imputed, FST-outlier-subset VCF -> PLINK bed for FEEMS.

# Usage: bash 5_vcf_to_plink.sh output/fst_splits/rep01/LOSH_topfst.vcf.gz input/feems/LOSH.feems_meta.csv input/feems/LOSH.top_pairwiseFST.second

set -euo pipefail
VCF=$1        # imputed VCF with breeding + nonbreeding birds
META=$2       # sample_id,group,lon,lat
OUT=$3        # output prefix

# samples in meta (skip header)
tail -n +2 "$META" | cut -d, -f1 > ${OUT}.keep.txt

# biallelic SNPs, only the birds in meta, unique IDs, GT hard calls
bcftools view -S ${OUT}.keep.txt --force-samples -m2 -M2 -v snps "$VCF" -Ou \
 | bcftools annotate --set-id '%CHROM:%POS' -Oz -o ${OUT}.tmp.vcf.gz
bcftools index -f ${OUT}.tmp.vcf.gz

# --double-id so IID == VCF sample name (what the python script matches on)
plink --vcf ${OUT}.tmp.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
      --make-bed --out ${OUT}

# Report residual missingness (should be ~0 after imputation)
plink --bfile ${OUT} --allow-extra-chr --missing --out ${OUT}.miss >/dev/null
awk 'NR>1 && $6>0' ${OUT}.miss.imiss | wc -l | xargs echo "samples with any missing GT:"

# OPTIONAL: outlier SNPs often sit in LD blocks; thinning stops a few regions
# dominating the likelihood. Compare assignments with/without.
 plink --bfile ${OUT} --allow-extra-chr --indep-pairwise 50 5 0.5 --out ${OUT}.ld
 plink --bfile ${OUT} --allow-extra-chr --extract ${OUT}.ld.prune.in --make-bed --out ${OUT}.exclLD

rm -f ${OUT}.tmp.vcf.gz*
echo "wrote ${OUT}.bed/.bim/.fam"

