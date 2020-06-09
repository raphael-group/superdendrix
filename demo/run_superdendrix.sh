#!/bin/bash

source activate superdendrix-env
cd ../

# downlaod data
cd data
snakemake all
cd ../

# define parameters
RANDSEED=2019
DIFF_DEP="data/profiles/CERES_scores_2C.tsv"
FEATURES="data/features/Q1_CERES_full_MUT.tsv"
GMM_OUTPUT="data/profiles/Q1_GMM"
CYCLE_FEATURES=10
OUTPUT_RANDMAT="data/features/rand_mat/"
PRE=0

THREADS=4
DIFF_DEP_SCORES="data/profiles/Q1_GMM.tsv"
GENE="KRAS"
FEATURES_ONCOKB="data/features/Q1_CERES_oncoKB_MUT.tsv"
CYCLE_SD=10
CYCLE_CP=10
DIRECTION="negative"
SETSIZE=3
OUTPUT_SD="data/results/KRAS.txt"

# compute CERES z-scores and identify six-sigma genes

"check packages - revealerpy-env"
echo "1"
Rscript utils/compute_CERES_zscores.R

# identify differential dependencies among the six-sigma genes
echo "2"
Rscript src/fit_tmm.R

# score the differential dependencies
echo "3"
"start from here"
python utils/fit_gmm.py -pf ${DIFF_DEP} -o ${GMM_OUTPUT}

# maf annotation with OncoKB
echo "4"
Rscript src/oncokb_maf_annotator.R


# generate randomized feature matrices
python utils/generate_null_matrices.py -m ${FEATURES} -p ${CYCLE_FEATURES} -o ${OUTPUT_RANDMAT} -pre ${PREFIX} -rs ${RANDSEED}

# identify the optimal association and perform statistical evaluation
python src/superdendrix.py -t ${THREADS} -T ${DIFF_DEP_SCORES} -Tc ${GENE} -m ${FEATURES} -p ${CYCLE_SD} -cp ${CYCLE_CP} -d ${DIRECTION} -k ${SETSIZE} -nm ${OUTPUT_RANDMAT} -rs ${RANDSEED} -x -curve -o ${OUTPUT_SD}
