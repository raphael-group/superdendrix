#!/usr/bin/env bash

# CONSTRUCT DIRECTORIES
RAW_DIR="data/raw/"
FEATURE_DIR="data/features/"
PROFILE_DIR="data/profiles/"

mkdir -p $FEATURE_DIR
mkdir -p $PROFILE_DIR

#################################################################
# DIFFERENTIAL DEPENDENCY COMPUTATION
#################################################################
## Normalize dependency scores into z-scores
DEPENDENCY_SCORES="$RAW_DIR"toy_dependency_scores.csv
ZSCORES="$PROFILE_DIR"dependency_zscores.csv
Rscript utils/verified/compute_zscores.R $DEPENDENCY_SCORES $ZSCORES

# #################################################################
# ## Fit t-mixture distributions to the z-scores
# INPUT_SCORES="$ZSCORES"
# FIT_RESULTS="$PROFILE_DIR"tmm_summary.tsv
# Rscript utils/verified/fit_tmm_20Q2.R $INPUT_SCORES $FIT_RESULTS

# #################################################################
# ## Fit gaussian mixture distributions to the z-scores
# INPUT_SCORES="$ZSCORES"
# GMM_DIR="$PROFILE_DIR"GMM/
# mkdir -p $GMM_DIR
# python utils/verified/fit_gmm.py -pf ${ZSCORES} -o ${GMM_DIR}

# ##################################################################
# ## Ensure monotonicity between z-scores and 2C scores
# GMM_SUMMARY="$GMM_DIR"gmm_profile_summary.tsv
# TWOC_SCORES="$GMM_DIR"twoC_scores.tsv
# DEPENDENCY_SCORES="$ZSCORES"
# OUTPUT_DIRECTORY="$PROFILE_DIR"

# python utils/verified/monotonicity.py -gmm $GMM_SUMMARY -twoc $TWOC_SCORES \
#     -raw $DEPENDENCY_SCORES \
#     -o $OUTPUT_DIRECTORY
