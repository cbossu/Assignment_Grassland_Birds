# Assignment_Grassland_Birds
Code to impute genotypes from an ANGSD beagle files, compute pairwise FST between groups and use FEEMS to generate a migration surface and assign nonbreeding birds to latitude/longitude coordinates

The first code, 1_beagle3_impute.sbatch uses beagle 3.1 to impute genotype probabilities from ANGSD into a vcf. 
The second code takes results from one to create a vcf and then to run pairwise FST between groups, take top FST outliers, unliked loci to create a -bfile in plink. If you already have a vcf file, you can skip the first script and first part of second script to create the -bfile from plink to run feems.