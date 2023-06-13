#!/bin/sh


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use Anaconda3

source activate ldsc


gwis_dir=../data/processed/gwis
ldref_dir=../data/processed/ld_ref
ldsc_datadir=../data/processed/ldsc


bm1=$1
bm2=$2


working_dir=${bm1}_${bm2}_ldsc_dir


# Munge summary statistics to prep for LDSC
mkdir -p ${working_dir}
for bm in ${bm1} ${bm2}; do 
cut  -f 2- ${gwis_dir}/${bm}_MD_gwis_merged > ${working_dir}/${bm}_ss
../opt/ldsc/munge_sumstats.py \
        --sumstats ${working_dir}/${bm}_ss \
        --merge-alleles ../data/raw/ldsc/w_hm3.snplist \
        --snp RSID \
        --N-col N_Samples \
        --a1 Effect_Allele \
        --a2 Non_Effect_Allele \
        --p P_Value_Marginal \
        --frq AF \
        --signed-sumstats Beta_Marginal,0 \
        --chunksize 500000 \
        --out ${working_dir}/${bm}
rm ${working_dir}/${bm}_ss
done


# Run LDSC to calculate genetic correlations
../opt/ldsc/ldsc.py \
	--rg ${working_dir}/${bm1}.sumstats.gz,${working_dir}/${bm2}.sumstats.gz \
	--ref-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--w-ld ${ldref_dir}/ukb_20k_hg19_withCM \
	--out ${ldsc_datadir}/${bm1}_${bm2}

rm -r ${working_dir}

conda deactivate
