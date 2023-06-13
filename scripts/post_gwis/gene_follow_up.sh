#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=2:00:00

#$ -cwd
#$ -j y


chr=$1
range=$2
gene_symbol=$3
biomarkers=$4
model=$5

geno_dir=../data/processed/genotypes
mkdir -p ${geno_dir}


source /broad/software/scripts/useuse
use R-4.0
use GCC-5.2


# Subset genome-wide summary stats to the gene region
for bm in ${biomarkers}; do
R --vanilla <<EOF
library(tidyverse)
pos_vec <- as.numeric(str_split(gsub(".*:", "", "${range}"), "-", simplify=TRUE))
read_tsv("../data/processed/gwis/${bm}_${model}_gwis_chr${chr}") %>%
  filter(POS >= pos_vec[1], POS <= pos_vec[2]) %>%
  write_tsv("../data/processed/gwis/${bm}_${model}_gwis_${gene_symbol}_subset")
EOF
done

# Use qctool v2 to subset (by position) to variants in the gene
qctool="~/kw/opt/qctool_v2.0.6-CentOS\ Linux7.3.1611-x86_64/qctool"
eval "${qctool}" \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
	-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
	-incl-range "${range}" \
	-og ${geno_dir}/${gene_symbol}_genotypes.bgen \
	-os ${geno_dir}/${gene_symbol}_genotypes.sample

# Convert genotype file to PLINK2 format
plink2=~/kw/opt/plink2
${plink2} \
	--bgen ${geno_dir}/${gene_symbol}_genotypes.bgen ref-first \
	--sample ${geno_dir}/${gene_symbol}_genotypes.sample \
	--make-pgen \
	--out ${geno_dir}/${gene_symbol}_genotypes
${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--export vcf bgz vcf-dosage=DS-force id-paste=iid \
	--out ${geno_dir}/${gene_symbol}_genotypes 
