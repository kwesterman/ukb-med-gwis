#!/bin/bash


biomarker=$1

source_datadir=../data/raw/westerman_2022
target_datadir=../data/processed/prior_evidence
annovardir=../opt/annovar

me_inputfile=${source_datadir}/main_effect_meta_analysis/${biomarker}_MA_1.tbl_nom
vqtl_inputfile=${source_datadir}/vqtl_meta_analysis/${biomarker}_MA_1.tbl_nom
annofile=${source_datadir}/all_chr.bim

# Create ANNOVAR input files

R --vanilla <<EOF
library(tidyverse)
bim_df <- read_tsv("${annofile}", col_names=c("CHR", "SNP", "CM", "POS", "ALT", "REF")) %>%
  select(SNP, CHR, POS)
read_tsv("${me_inputfile}") %>%
  inner_join(bim_df, by="SNP") %>%
  select(SNP, CHR, POS, REF, ALT, P) %>%
  select(CHR, start=POS, end=POS, REF, ALT, SNP, P) %>%
  write_tsv("${target_datadir}/${biomarker}_me.avinput")
read_tsv("${vqtl_inputfile}") %>%
  inner_join(bim_df, by="SNP") %>%
  select(SNP, CHR, POS, REF, ALT, P) %>%
  select(CHR, start=POS, end=POS, REF, ALT, SNP, P) %>%
  write_tsv("${target_datadir}/${biomarker}_vqtl.avinput")
EOF

# Run ANNOVAR

perl $annovardir/annotate_variation.pl \
        -out ${datadir}/${biomarker}_me \
        -build hg19 \
        ${target_datadir}/${biomarker}_me.avinput \
        ${annovardir}/humandb/

perl $annovardir/annotate_variation.pl \
        -out ${datadir}/${biomarker}_vqtl \
        -build hg19 \
        ${target_datadir}/${biomarker}_vqtl.avinput \
        ${annovardir}/humandb/
