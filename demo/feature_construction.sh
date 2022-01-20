#!/usr/bin/env bash

# CONSTRUCT DIRECTORIES
RAW_DIR="data/raw/"
FEATURE_DIR="data/features/"
PROFILE_DIR="data/profiles/"

mkdir -p $FEATURE_DIR
mkdir -p $PROFILE_DIR

#################################################################
# FEATURE CONSTRUCTION
#################################################################
## ONCOKB ANNOTATION
# echo "Start annotating the MAF file with OncoKB database"

# TOKEN="12cd9caf-9fb8-4098-80fd-9de027fe6b9a" # OncoKB API Token

MAF_INPUT="$RAW_DIR"CCLE_mutations_avana.tsv # Input mutations file
MAF_ONCOKB="$RAW_DIR"CCLE_mutations_avana_oncoKB.tsv # Name of output file
# python oncokb-annotator/MafAnnotator.py -i $MAF_INPUT -o $MAF_ONCOKB \
#     -b $TOKEN -q HGVSp_Short

#################################################################
## Feature construction
echo "Construct OncoKB features"
INPUT_MAF="$MAF_ONCOKB" # Use the previously annotaed mutations as input
SAMPLE_INFO="$RAW_DIR"sample_info_avana.tsv
# Rscript utils/verified/oncokb_maf_annotator_20Q2.R $INPUT_MAF \
#     $SAMPLE_INFO $FEATURE_DIR

## Merge Inactivating and Other mutations in the same gene
INPUT_FEATURES="$FEATURE_DIR"features_all_MUT.tsv
MERGED_FEATURES="$FEATURE_DIR"features_all_MUT_IOmerged.tsv
python utils/merge_I_O.py $INPUT_FEATURES $MERGED_FEATURES

INPUT_FEATURES="$FEATURE_DIR"/features_oncokb_mut.tsv
MERGED_FEATURES="$FEATURE_DIR"/features_oncokb_mut_IOmerged.tsv
python utils/merge_I_O.py $INPUT_FEATURES $MERGED_FEATURES

# #################################################################
# ## Generate randomized matrices (we use all of the mutations 
# ## in each cell line to account for background mutation rate)
# NUM_MATRICES=10    # number of randomized matrices to construct
# FEATURES="$MERGED_FEATURES"
# RAND_MAT_DIR="$FEATURE_DIR"rand_mat/
# mkdir -p ${RAND_MAT_DIR}

# echo Generate "$NUM_MATRICES" randomized matrices
# python utils/generate_null_matrices.py -m $FEATURES -p $NUM_MATRICES \
#     -o $RAND_MAT_DIR
