#!/bin/bash

# demo for analysis of CERES dataset using SuperDendrix algorithm

# source activate {your-env}  # activate the conda environment with python and r dependencies
source activate superdendrix-env

# download data
cd ../data

snakemake all

# make directories
mkdir profiles features
mkdir features/rand_mat
mkdir results

cd ../

# define parameters
RANDSEED=2019
DIFF_DEP="data/profiles/CERES_scores_two-comp.tsv"
FEATURES="data/features/Q1_CERES_full_MUT.tsv"
GMM_OUTPUT="data/profiles/Q1_2C_scores"
CYCLE_FEATURES=10
OUTPUT_RANDMAT="data/features/rand_mat/"
PREFIX=0

THREADS=4
DIFF_DEP_SCORES="data/profiles/Q1_2C_scores.tsv"
GENE="ARID1B (57492)"
#GENE="ACOX1 (51)"
FEATURES_ONCOKB="data/features/Q1_CERES_oncoKB_MUT.tsv"
CYCLE_SD=10
CYCLE_CP=10
DIRECTION="negative"
SETSIZE=3
OUTPUT_SD="data/results/ARID1B.txt"
FEATURE_LIST="data/features/Q1_CERES_oncoKB_feature_list.csv"
# compute CERES z-scores and identify six-sigma genes
##Rscript utils/compute_CERES_zscores.R

# identify differential dependencies among the six-sigma genes
##Rscript utils/fit_tmm.R

# score the differential dependencies
##python utils/fit_gmm.py -pf ${DIFF_DEP} -o ${GMM_OUTPUT}

# maf annotation with OncoKB
Rscript utils/oncokb_maf_annotator.R


# generate randomized feature matrices
python utils/generate_null_matrices.py -m ${FEATURES} -p ${CYCLE_FEATURES} -o ${OUTPUT_RANDMAT} -pre ${PREFIX} -rs ${RANDSEED}

# identify the optimal association and perform statistical evaluation
python src/superdendrix.py -t ${THREADS} -T ${DIFF_DEP_SCORES} -Tc "${GENE}" -m ${FEATURES_ONCOKB} -p ${CYCLE_SD} -cp ${CYCLE_CP} -d ${DIRECTION} -gf ${FEATURE_LIST} -k ${SETSIZE} -nm ${OUTPUT_RANDMAT} -rs ${RANDSEED} -x -curve -o ${OUTPUT_SD}
