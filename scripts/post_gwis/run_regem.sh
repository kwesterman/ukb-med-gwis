#!/bin/sh


#$ -l h_vmem=30G
#$ -l h_rt=1:00:00

#$ -j y
#$ -cwd


sumstats_fn=$1
new_exp=$2
new_intCov=$3
new_tag=$4

REGEM_dir=../data/processed/regem
mkdir -p ${REGEM_dir}
new_fn=${REGEM_dir}/${new_tag}_derived


source /broad/software/scripts/useuse
use .boost-1.70.0
use .mkl-2019.3.199

REGEM=../opt/REGEM/src/REGEM

${REGEM} \
	--input-file ${sumstats_fn} \
	--exposure-names ${new_exp} \
	--int-covar-names ${new_intCov} \
	--output-style minimum \
	--out ${new_fn}
