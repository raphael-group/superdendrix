from os.path import join

################################################################################
# SETTINGS, FILES, AND DIRECTORIES
################################################################################
# Directories

DATA_DIR = 'data/'

RAW_DIR = 'raw'
FEATURES_DIR = 'features'
PROFILES_DIR = 'profiles'



# URLs
##### Example dataset ####
CCLE_MAF_URL = 'https://www.dropbox.com/s/bfvvm5iwhr348jy/mutations.tsv?dl=0'
CELL_LINE_INFO_URL = 'https://www.dropbox.com/s/h56qxjvxlqs6re6/sample_info.tsv?dl=0'
DEPENDENCY_SCORES_URL = 'https://www.dropbox.com/s/wew164e9zxkvw1f/dependency_scores.csv?dl=0'

# Files
CCLE_MAF = join(DATA_DIR+RAW_DIR, 'mutations.tsv')
DEPENDENCY_SCORES = join(DATA_DIR+RAW_DIR, 'dependency_scores.csv')
CELL_LINE_INFO = join(DATA_DIR+RAW_DIR, 'sample_info.tsv')

################################################################################
# COMMANDS
################################################################################

# General
rule all:
    input:
        CCLE_MAF,
        DEPENDENCY_SCORES,
        CELL_LINE_INFO,

#Download
rule download_ccle_maf:
    params:
        url=CCLE_MAF_URL
    output:
        CCLE_MAF
    shell:
        'wget -O {output} {params.url}'

rule download_dependency_scores:
    params:
        url=DEPENDENCY_SCORES_URL
    output:
        DEPENDENCY_SCORES
    shell:
        'wget -O {output} {params.url}'

rule download_ceres_cell_line_info:
    params:
        url=CELL_LINE_INFO_URL
    output:
        CELL_LINE_INFO
    shell:
        'wget -O {output} {params.url}'
