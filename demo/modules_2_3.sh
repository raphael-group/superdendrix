#!/usr/bin/env bash

# CONSTRUCT DIRECTORIES

#################################################################
# RUN SUPERDENDRIX ILP
#################################################################
# Identify associations between differential dependency and
# mutually exclusive alteration features

# Input data
FEATURES='data/features/features_oncokb_MUT_IOmerged.tsv'  #OncoKB features
SCORES='data/profiles/GMM/monotonic_twoC_scores.tsv'  # 2C dependency scores
RANDOM_MATRICES='data/features/rand_mat/' # randomized mutation matrices

# Input parameters
GENE='BRAF_(673)' # Name of the differential dependency. Must be a column in SCORES
DIRECTION="negative" # Direction of dependency. e.g, negative: increased dependency, positive: decreased dependency

K=3 # Maximum number of mutations in a set
PERM_MODEL=10 # Number of permutations for model selection
PERM_ASSOCIATION=10 # Number of permutations for feature association test



# Set output directory and filename
OUTPUT_DIR='data/results/'
OUTPUT_FILENAME="test.tsv"
mkdir -p $OUTPUT_DIR


python src/superdendrix.py -T $SCORES -Tc $GENE -d $DIRECTION \
    -m $FEATURES -k $K \
    -nm $RANDOM_MATRICES -cp $PERM_MODEL -p $PERM_ASSOCIATION \
    -o $OUTPUT_DIR$OUTPUT_FILENAME
