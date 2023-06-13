#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=48:00:00

#$ -cwd
#$ -j y


chr=$1
range=$2
gene_symbol=$3
geno_dir=../data/processed/genotypes
mkdir -p ${geno_dir}


# Use qctool v2 to subset by position
source /broad/software/scripts/useuse
use GCC-5.2
qctool="~/kw/opt/qctool_v2.0.6-CentOS\ Linux7.3.1611-x86_64/qctool"

eval "${qctool}" \
	-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
	-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
	-incl-range "${range}" \
	-og ${geno_dir}/${gene_symbol}_genotypes.bgen \
	-os ${geno_dir}/${gene_symbol}_genotypes.sample

# Convert to PLINK2 format
plink2=~/kw/opt/plink2
${plink2} \
	--bgen ${geno_dir}/${gene_symbol}_genotypes.bgen ref-first \
	--sample ${geno_dir}/${gene_symbol}_genotypes.sample \
	--make-pgen \
	--out ${geno_dir}/${gene_symbol}_genotypes

${plink2} \
	--pfile ${geno_dir}/${gene_symbol}_genotypes \
	--export vcf bgz vcf-dosage=DS-force \
	--out ${geno_dir}/${gene_symbol}_genotypes 
