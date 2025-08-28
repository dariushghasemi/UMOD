#!/bin/bash

#SBATCH --job-name dosage
#SBATCH --output %j_extract_dosage.log
#SBATCH --partition cpuq
#SBATCH --cpus-per-task 1
#SBATCH --mem 13G
#SBATCH --time 00-01:00:00


# path to genotype
genotype=/scratch/dariush.ghasemi/projects/chr16.dose.vcf.gz
# initial VCF was: /scratch/ekoenig/CHRIS_CORRECTED/10K/Imputed/TOPMedR2/20210409/chr16.vcf.gz

# variants file
snps=/scratch/dariush.ghasemi/projects/UMOD/data/umod_hg38.snps
# coordinates of a geneomic region (also position of interested SNP)
# -r: region can also be defined as chr16_20381010-20381010

# path to output
ofile=/scratch/dariush.ghasemi/projects/UMOD/data/umod_hg38.dosage


# load libraries on server
source /exchange/healthds/singularity_functions

# Extracting dosage level of UMOD lead SNP in CHRIS
bcftools query -f '[%SAMPLE\t%ID\t%DS\n]' $genotype -R $snps -o $ofile

# for full characteristics of the variant, use below term for the query
#-f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%GT\t%DS\t%HDS\t%GP\n]'
