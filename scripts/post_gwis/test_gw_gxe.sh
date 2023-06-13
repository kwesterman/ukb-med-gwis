#!/bin/sh


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use Anaconda3



gwis_dir=../data/processed/gwis
ldref_dir=../data/processed/ld_ref
gw_gxe_dir=../data/processed/gw_gxe


bm=$1
exp=$2


working_dir=${bm}_${exp}_ldsc_dir


source activate ldsc

# Munge summary statistics to prep for LDSC
mkdir -p ${working_dir}
cut  -f 2- ${gwis_dir}/${bm}_${exp}_gwis_merged > ${working_dir}/${bm}_${exp}_ss
../opt/ldsc/munge_sumstats.py \
        --sumstats ${working_dir}/${bm}_${exp}_ss \
        --merge-alleles ../data/raw/ldsc/w_hm3.snplist \
        --snp RSID \
        --N-col N_Samples \
        --a1 Effect_Allele \
        --a2 Non_Effect_Allele \
        --p P_Value_Interaction \
        --frq AF \
        --signed-sumstats Beta_G-${exp},0 \
        --chunksize 500000 \
        --out ${working_dir}/${bm}_${exp}
rm ${working_dir}/${bm}_${exp}_ss

# Run LDSC to generate genome-wide GxE heritability estimate
../opt/ldsc/ldsc.py \
	--h2 ${working_dir}/${bm}_${exp}.sumstats.gz \
	--ref-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--w-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--out ${gw_gxe_dir}/${bm}_${exp}_ldsc_h2

conda deactivate

source activate pigeon

# Run PIGEON to generate genome-wide GxE heritability estimate
../opt/PIGEON/pigeon.py \
	--gxe-var \
	--gxe-sumstats ${working_dir}/${bm}_${exp}.sumstats.gz \
	--ref-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--w-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--out ${gw_gxe_dir}/${bm}_${exp}_pigeon_h2

conda deactivate

rm -r ${working_dir}
