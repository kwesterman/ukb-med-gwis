#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=5G
#$ -l h_rt=36:00:00

#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -j y
#$ -cwd


pheno=$1

chr=$SGE_TASK_ID


gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5)
med_scoreBygPC_arr=( "${gPC_arr[@]/#/med_scoreBy}" )
minimal_covars="age sex ageBysex ${gPC_arr[@]} ${minimal_covars} ${med_scoreBygPC_arr[@]}"

md_components="VEG LEGUMES FRUIT NUTS FISH WHGRAIN MUFA2SFA REDPRMEAT ALC"


singularity_dir=~/kw/singularity
singularity exec \
	-B ../data/processed:/data \
	-B /broad/ukbb/imputed_v3:/bgendir \
	-B /humgen/florezlab/UKBB_app27892:/sampledir \
	${singularity_dir}/gem-v1.3-workflow.simg \
	/bin/bash <<EOF

/GEM/GEM \
	--bgen /bgendir/ukb_imp_chr${chr}_v3.bgen \
	--sample /sampledir/ukb27892_imp_chrAUT_v3_s487395.sample \
	--pheno-file /data/ukb_phenos_unrelated.csv \
	--sampleid-name id \
	--pheno-name ${pheno} \
	--pheno-type 0 \
	--exposure-names ${md_components} \
	--covar-names ${minimal_covars} \
	--delim , \
	--maf 0.01 \
	--missing-value NA \
	--robust 1 \
	--threads 8 \
	--output-style minimum \
	--out /data/gwis/${pheno}_MDcomponents_gwis_chr${chr}

EOF
