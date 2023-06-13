#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=00:30:00
#$ -j y
#$ -cwd


pheno=$1

geno_dir=$(readlink -f ../data/processed/genotypes)
res_dir=$(readlink -f ../data/processed/ukb_gxe_res)
mkdir -p ${res_dir}

gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5)
med_scoreBygPC_arr=( "${gPC_arr[@]/#/med_scoreBy}" )
minimal_covars="age sex ageBysex ${gPC_arr[@]} ${minimal_covars} ${med_scoreBygPC_arr[@]}"
#minimal_covars="age sex ageBysex ${gPC_arr[@]}"

echo "" > ${res_dir}/${pheno}_minimal

for chr in {1..22} X; do

singularity_dir=~/kw/singularity
singularity exec \
	-B ../data/processed:/data \
	-B ${geno_dir}:/geno_dir \
	-B ${res_dir}:/res_dir \
	${singularity_dir}/gem-v1.3-workflow.simg \
	/bin/bash <<EOF

/GEM/GEM \
	--bgen /geno_dir/target_genotypes_chr${chr}.bgen \
	--sample /geno_dir/target_genotypes_chr${chr}.sample \
	--pheno-file /data/ukb_phenos_unrelated.csv \
	--sampleid-name id \
	--pheno-name ${pheno} \
	--pheno-type 0 \
	--exposure-names med_score \
	--covar-names ${minimal_covars} \
	--delim , \
	--maf 0.0001 \
	--missing-value NA \
	--robust 1 \
	--threads 1 \
	--output-style minimum \
	--out /res_dir/${pheno}_minimal_chr${chr}

if [[ "$chr" -eq 1 ]]; then
	head -1 /res_dir/${pheno}_minimal_chr${chr} >> /res_dir/${pheno}_minimal
fi

tail -n +2 /res_dir/${pheno}_minimal_chr${chr} >> /res_dir/${pheno}_minimal

EOF

done

