#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=10G
#$ -l h_rt=00:30:00
#$ -j y
#$ -cwd


pheno=$1
exposure=$2

geno_dir=$(readlink -f ../data/processed/genotypes)
res_dir=$(readlink -f ../data/processed/ukb_gxe_res)
mkdir -p ${res_dir}

gPC_arr=(gPC1 gPC2 gPC3 gPC4 gPC5)
med_scoreBygPC_arr=( "${gPC_arr[@]/#/med_scoreBy}" )
minimal_covars="age sex ageBysex ${gPC_arr[@]} ${minimal_covars} ${med_scoreBygPC_arr[@]}"

echo "" > ${res_dir}/${pheno}_minimal

singularity_dir=~/kw/singularity
singularity exec \
        -B ../data/processed:/data \
        -B ${geno_dir}:/geno_dir \
        -B ${res_dir}:/res_dir \
        ${singularity_dir}/gem-v1.3-workflow.simg \
        /bin/bash <<EOF

/GEM/GEM \
        --bgen /geno_dir/slc_genotypes.bgen \
        --sample /geno_dir/slc_genotypes.sample \
        --include-snp-file /geno_dir/followup_rsids.txt \
        --pheno-file /data/ukb_phenos_unrelated.csv \
        --sampleid-name id \
        --pheno-name ${pheno} \
        --pheno-type 0 \
        --exposure-names ${exposure} \
        --covar-names ${minimal_covars} \
        --delim , \
        --maf 0.0001 \
        --missing-value NA \
        --robust 1 \
        --threads 1 \
        --output-style minimum \
        --out /res_dir/${pheno}_minimal_followup_exp_${exposure}

EOF
