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

# path to output
ofile=/scratch/dariush.ghasemi/projects/UMOD/data/chr16_20381010_rs77924615.txt

# coordinates of a geneomic region (also position of interested SNP)
region=chr16:20381010   # hg38 for UMOD hit in CHRIS (rs77924615)
#region can also be defined as chr16_20381010-20381010

# load libraries on server
source /exchange/healthds/singularity_functions

# Extracting dosage level of UMOD lead SNP in CHRIS
bcftools query -f '[%SAMPLE\t%ID\t%DS\n]' $genotype -r $region -o $ofile

# for full characteristics of the variant, use below term for the query
#-f '[%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%GT\t%DS\t%HDS\t%GP\n]'
