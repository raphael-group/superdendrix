#!/usr/bin/env bash

# CONSTRUCT DIRECTORIES
RAW_DIR="data/raw/" # directory for original dataset

FEATURE_DIR="data/features/" # directory for oncokb features
PROFILE_DIR="data/profiles/" # directory for differential dependencies

mkdir -p $FEATURE_DIR
mkdir -p $PROFILE_DIR

#################################################################
# DIFFERENTIAL DEPENDENCY COMPUTATION
#################################################################
## Normalize dependency scores into z-scores
echo "Normalize dependency scores into z-scores"
DEPENDENCY_SCORES="$RAW_DIR"dependency_scores.csv
ZSCORES="$PROFILE_DIR"dependency_zscores.tsv

Rscript src/compute_zscores.R $DEPENDENCY_SCORES $ZSCORES

#################################################################
## Fit t-mixture distributions to the z-scores
echo "Fitting t-mixture model to z-scores"
INPUT_SCORES="$ZSCORES"
FIT_RESULTS="$PROFILE_DIR"tmm_summary.tsv

Rscript src/fit_tmm.R $INPUT_SCORES $FIT_RESULTS

#################################################################
## Fit gaussian mixture distributions to the z-scores
echo "fitting Gaussian mixture model to z-scores"
INPUT_SCORES="$ZSCORES"
GMM_DIR="$PROFILE_DIR"GMM/
mkdir -p $GMM_DIR

python src/fit_gmm.py -pf ${ZSCORES} -o ${GMM_DIR}

##################################################################
## Ensure monotonicity between z-scores and 2C scores
echo "Correct the scores for monotonicity"
GMM_SUMMARY="$GMM_DIR"gmm_profile_summary.tsv
TWOC_SCORES="$GMM_DIR"twoC_scores.tsv
DEPENDENCY_SCORES="$ZSCORES"
OUTPUT_DIRECTORY="$GMM_DIR"

python src/monotonicity.py -gmm $GMM_SUMMARY -twoc $TWOC_SCORES \
    -raw $DEPENDENCY_SCORES \
    -o $OUTPUT_DIRECTORY


#################################################################
# FEATURE CONSTRUCTION
#################################################################
## ONCOKB ANNOTATION
echo "Start annotating the MAF file with OncoKB database"

TOKEN="" # Obtain and add your OncoKB API Token from https://www.oncokb.org/

MAF_INPUT="$RAW_DIR"mutations.tsv # Input mutations file
MAF_ONCOKB="$RAW_DIR"mutations_oncoKB.tsv # Name of output file
python oncokb-annotator/MafAnnotator.py -i $MAF_INPUT -o $MAF_ONCOKB \
    -b $TOKEN -q HGVSp_Short

#################################################################
## Feature construction
echo "Construct OncoKB features"
INPUT_MAF="$MAF_ONCOKB" # Use the previously annotaed mutations as input
SAMPLE_INFO="$RAW_DIR"sample_info.tsv
Rscript src/construct_oncokb_features.R $INPUT_MAF \
    $SAMPLE_INFO $FEATURE_DIR

## Merge Inactivating and Other mutations in the same gene
INPUT_FEATURES="$FEATURE_DIR"features_all_MUT.tsv
MERGED_FEATURES="$FEATURE_DIR"features_all_MUT_IOmerged.tsv
python src/merge_I_O.py -i $INPUT_FEATURES -o $MERGED_FEATURES

INPUT_FEATURES="$FEATURE_DIR"features_oncokb_MUT.tsv
MERGED_FEATURES="$FEATURE_DIR"features_oncokb_MUT_IOmerged.tsv
python src/merge_I_O.py -i $INPUT_FEATURES -o $MERGED_FEATURES

#################################################################
## Generate randomized matrices (we use all of the mutations 
## in each cell line to account for background mutation rate)
NUM_MATRICES=10    # number of randomized matrices to construct
FEATURES="$MERGED_FEATURES"
RAND_MAT_DIR="$FEATURE_DIR"rand_mat/
mkdir -p ${RAND_MAT_DIR}

echo Generate "$NUM_MATRICES" randomized matrices
python src/generate_null_matrices.py -m $FEATURES -p $NUM_MATRICES \
    -o $RAND_MAT_DIR
