#!/bin/bash


#$ -l h_vmem=20G
#$ -l h_rt=48:00:00

#$ -cwd
#$ -j y


variant_file=../data/processed/target_gene_variant_list.tsv
geno_dir=../data/processed/genotypes
mkdir -p ${geno_dir}


# Generate file of variant targets
cut -f 3 ${variant_file} > ${geno_dir}/target_rsids.txt


# Use qctool v2 to subset 
source /broad/software/scripts/useuse
use GCC-5.2
qctool="../../opt/qctool_v2.0.6-CentOS\ Linux7.3.1611-x86_64/qctool"

for chr in {1..22}; do
	eval "${qctool}" \
		-g /broad/ukbb/imputed_v3/ukb_imp_chr${chr}_v3.bgen \
		-s /humgen/florezlab/UKBB_app27892/ukb27892_imp_chrAUT_v3_s487395.sample \
		-incl-rsids ${geno_dir}/target_rsids.txt \
		-og ${geno_dir}/target_genotypes_chr${chr}.bgen \
		-os ${geno_dir}/target_genotypes_chr${chr}.sample
done

# Merge subsets
cp ${geno_dir}/target_genotypes_chr22.bgen ${geno_dir}/target_genotypes_base.bgen  # Initialize "base" .bgen file as chr22 .bgen
cp ${geno_dir}/target_genotypes_chr22.sample ${geno_dir}/target_genotypes_base.sample
for chr in {21..1}; do
	echo "Merging chromosome ${chr}..."
	eval "${qctool}" \
		-g ${geno_dir}/target_genotypes_base.bgen \
		-s ${geno_dir}/target_genotypes_base.sample \
		-merge-in ${geno_dir}/target_genotypes_chr${chr}.bgen ${geno_dir}/target_genotypes_chr${chr}.sample \
		-og ${geno_dir}/target_genotypes.bgen \
		-os ${geno_dir}/target_genotypes.sample
	cp ${geno_dir}/target_genotypes.bgen ${geno_dir}/target_genotypes_base.bgen  # Merged file becomes the new base for next round
	cp ${geno_dir}/target_genotypes.sample ${geno_dir}/target_genotypes_base.sample
done
rm ${geno_dir}/target_genotypes_base.bgen ${geno_dir}/target_genotypes_base.sample

# Convert to PLINK2 format
plink2=../../opt/plink2
${plink2} \
	--bgen ${geno_dir}/target_genotypes.bgen ref-first \
	--sample ${geno_dir}/target_genotypes.sample \
	--make-pgen \
	--out ${geno_dir}/target_genotypes
