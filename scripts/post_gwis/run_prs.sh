#!/bin/sh


#$ -l h_vmem=10G
#$ -l h_rt=1:00:00
#$ -cwd
#$ -j y


gwis_dir=../data/processed/gwis
ldref_datadir=../data/processed/ld_ref
prs_datadir=../data/processed/prs


bm1=$1
bm2=$2


working_dir=${bm1}_${bm2}_prs_dir
mkdir -p ${working_dir}


source /broad/software/scripts/useuse
use R-4.0


# Remove duplicates from summary stats
#R --no-save <<EOF
#library(tidyverse)
#read_tsv("${gwis_dir}/${bm1}_MD_gwis_merged") %>%
#  distinct(RSID, .keep_all=TRUE) %>%
#  write_tsv("${working_dir}/sumstats.tsv")
#EOF
awk -F'\t' '!seen[$2]++' ${gwis_dir}/${bm1}_MD_gwis_merged > ${working_dir}/sumstats.tsv

plink2=../opt/plink2.0/plink2
${plink2} \
	--bfile ${ldref_datadir}/ukb_20k_hg19_withCM \
	--rm-dup retain-mismatch list \
	--out ${working_dir}/dups

# Generate a P&T PRS using BM1 main effect summary statistics
plink19=../opt/plink1.9/plink
${plink19} \
	--bfile ${ldref_datadir}/ukb_20k_hg19_withCM \
	--exclude ${working_dir}/dups.rmdup.list \
	--clump ${working_dir}/sumstats.tsv \
	--clump-snp-field RSID \
	--clump-field P_Value_Marginal \
	--out ${working_dir}/prs

awk 'NR!=1 {print $3}' ${working_dir}/prs.clumped > ${working_dir}/prs.snps

# Calculate the resulting PRS in the same sample set
${plink19} \
	--bfile ${ldref_datadir}/ukb_20k_hg19_withCM \
	--score ${working_dir}/sumstats.tsv 2 6 9 header \
	--extract ${working_dir}/prs.snps \
	--out ${working_dir}/prs

# Test the association between the BM1 PRS and BM2
R --no-save <<EOF
library(tidyverse)
prs_df <- read_table("${working_dir}/prs.profile") %>%
  select(id=IID, score=SCORE)
phenos <- read_csv("../data/processed/ukb_phenos_unrelated.csv") %>%
  inner_join(prs_df, by="id")
lm_fit <- lm(${bm2} ~ score, data=phenos)
saveRDS(lm_fit, "${prs_datadir}/${bm1}_prs_${bm2}_lm_fit.rds")
EOF

rm -r ${working_dir}
