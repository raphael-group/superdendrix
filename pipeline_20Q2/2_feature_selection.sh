#!/usr/bin/env bash

# CONSTRUCT DIRECTORIES
RAW_DIR="data/raw/"
FEATURE_DIR="data/features/"
PROFILE_DIR="data/profiles/"

#################################################################
# FEATURE CONSTRUCTION
#################################################################
## ONCOKB ANNOTATION
echo "Start annotating the MAF file with OncoKB database"

TOKEN="12cd9caf-9fb8-4098-80fd-9de027fe6b9a" # OncoKB API Token

MAF_INPUT="$RAW_DIR"CCLE_mutations_avana.tsv # Input mutations file
MAF_ONCOKB="$RAW_DIR"CCLE_mutations_avana_oncoKB.tsv # Name of output file
python oncokb-annotator/MafAnnotator.py -i $MAF_INPUT -o $MAF_OUTPUT \
    -b $TOKEN -q HGVSp_Short

#################################################################
## Feature construction
echo "Construct OncoKB features"
INPUT_MAF=$"MAF_ONCOKB" # Use the previously annotaed mutations as input
SAMPLE_INFO="$RAW_DIR"sample_info_avana.tsv
Rscript utils/verified/oncokb_maf_annotator_20Q2.R $INPUT_MAF \
    $SAMPLE_INFO $FEATURE_DIR

## Merge Inactivating and Other mutations in the same gene
INPUT_FEATURES="$FEATURE_DIR"features_all_MUT.tsv
MERGED_FEATURES="$FEATURE_DIR"features_all_MUT_IOmerged.tsv
python utils/merge_I_O.py $IMAF $OMAF

INPUT_FEATURES="$FEATURE_DIR"/features_oncokb_mut.tsv
MERGED_FEATURES="$FEATURE_DIR"/features_oncokb_mut_IOmerged.tsv
python utils/merge_I_O.py $IMAF $OMAF

#################################################################
## Generate randomized matrices (we use all of the mutations 
## in each cell line to account for background mutation rate)

NUM_MATRICES=10          # number of randomized matrices to construct
FEATURES="$MERGED_FEATURES"
RAND_MAT_DIR="$FEATURE_DIR"rand_mat/
mkdir -p ${RAND_MAT_DIR}

python utils/generate_null_matrices.py -m $FEATURES -p $NUM_MATRICES \
    -o $RAND_MAT_DIR

#################################################################
# DIFFERENTIAL DEPENDENCY COMPUTATION
#################################################################
## Normalize dependency scores into z-scores
DEPENDENCY_SCORES="$RAW_DIR"toy_dependency_scores.csv
ZSCORES="$PROFILE_DIR"dependency_zscores.csv
Rscript utils/verified/compute_zscores.R $DEPENDENCY_SCORES $ZSCORES

#################################################################
## Fit t-mixture distributions to the z-scores
INPUT_SCORES="$ZSCORES"
FIT_RESULTS="$PROFILE_DIR"tmm_summary.tsv
Rscript utils/verified/fit_tmm_20Q2.R $INPUT_SCORES $FIT_RESULTS

#################################################################
## Fit gaussian mixture distributions to the z-scores
INPUT_SCORES="$ZSCORES"
GMM_DIR="$PROFILE_DIR"GMM/
mkdir -p $GMM_DIR
python utils/verified/fit_gmm.py -pf ${ZSCORES} -o ${GMM_DIR}

##################################################################
## Ensure monotonicity between z-scores and 2C scores
GMM_SUMMARY="$GMM_DIR"gmm_profile_summary.tsv
TWOC_SCORES="$GMM_DIR"twoC_scores.tsv
DEPENDENCY_SCORES="$ZSCORES"
OUTPUT_DIRECTORY="$PROFILE_DIR"

python utils/verified/monotonicity.py -gmm $GMM_SUMMARY -twoc $TWOC_SCORES \
    -raw $DEPENDENCY_SCORES \
    -o $OUTPUT_DIRECTORY

#################################################################
# RUN SUPERDENDRIX ILP
#################################################################
# Identify associations between differential dependency and
# mutually exclusive alteration features

ALT='MUT'
K=3 # Maximum number of mutations
PERM_MODEL=10000 # Number of permutations for model selection
PERM_ASSOCIATION=10000 # Number of permutations for feature association test


FEATURES='data/21Q1/features/IOmerged/CERES_OncoKB_'$ALT'_merged.tsv'	#Ceres alterations
SCORES='data/21Q1/profiles/GMM/corrected_CERES_zscores_2C_with_resp_2Cscores.tsv'
RANDOM_MATRICES='data/21Q1/features/rand_mat/'$ALT'/' #full mut NONSILENT

OUTPUT_DIR='experiments/21Q1/gmm/scrap/'
OUTPUT_FILENAME="test.tsv"
mkdir -p $OUTPUT_DIR

GENE='BRAF_(673)'
DIRECTION="negative"
python src/superdendrix.py -T $SCORES -Tc $GENE -d $DIRECTION \
    -m $FEATURES -k $K \
    -nm $NULLMAT -cp $PERM_MODEL -p $PERM_ASSOCIATION \
    -o $OUTPUT_DIR$OUTPUT_FILENAME