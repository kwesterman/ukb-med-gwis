#!/bin/bash


#$ -l h_rt=2:00:00

#$ -cwd
#$ -j y


rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-108/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz .
zcat 1000GENOMES-phase_3.vcf.gz | tail -n +67 | cut -f 1-5 > snp_annot_hg19.txt
#paste snp_annot_hg19.txt <(awk 'BEGIN {OFS=":"} {print $1,$2}' snp_annot_hg19.txt) > snp_annot_hg19_withLoci.txt
#cat snp_annot_hg19_withLoci.txt | sort -u -k3,3 | sort -u -k6,6 -V > snp_annot_hg19_nodups.txt  # Unique rsID and chr:pos locus
#rm snp_annot_hg19.txt snp_annot_hg19_withLoci.txt
