from os.path import join

################################################################################
# SETTINGS, FILES, AND DIRECTORIES
################################################################################
# Directories

DIR_20Q2 = 'data/20Q2/'

RAW_DIR = 'raw'
FEATURES_DIR = 'features'
#FEATURES_RM_DIR = 'features/rand_mat'
PROFILES_DIR = 'profiles'
RESULTS_DIR = 'results'



# URLs
##### 20Q2 ####
CCLE_MAF_URL = 'https://ndownloader.figshare.com/files/22629110'
CERES_SCORES_URL = 'https://ndownloader.figshare.com/files/22629068'
CELL_LINE_INFO_URL = 'https://ndownloader.figshare.com/files/22629137'

# Files
CCLE_MAF = join(DIR_20Q2+RAW_DIR, 'mutations.csv')
CERES_SCORES = join(DIR_20Q2+RAW_DIR, 'gene_effect.csv')
CELL_LINE_INFO = join(DIR_20Q2+RAW_DIR, 'sample_info.csv')

################################################################################
# COMMANDS
################################################################################

# General
rule all:
    input:
        CCLE_MAF,
        CERES_SCORES,
        CELL_LINE_INFO,

#Download
rule download_ccle_maf:
    params:
        url=CCLE_MAF_URL
    output:
        CCLE_MAF
    shell:
        'wget -O {output} {params.url}'

rule download_ceres_scores:
    params:
        url=CERES_SCORES_URL
    output:
        CERES_SCORES
    shell:
        'wget -O {output} {params.url}'

rule download_ceres_cell_line_info:
    params:
        url=CELL_LINE_INFO_URL
    output:
        CELL_LINE_INFO
    shell:
        'wget -O {output} {params.url}'