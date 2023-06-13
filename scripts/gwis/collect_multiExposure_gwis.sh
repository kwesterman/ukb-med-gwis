#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=00:30:00

#$ -cwd
#$ -j y


fn_prefix=$1


head -2 ${fn_prefix}_chr1 > ${fn_prefix}_merged
for chr in {1..22}; do
	echo "${fn_prefix}_chr${chr}..."
	tail -n +3 ${fn_prefix}_chr${chr} >> ${fn_prefix}_merged
done
