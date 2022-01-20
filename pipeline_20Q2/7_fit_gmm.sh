#!/bin/bash

#SBATCH --partition raphael
#SBATCH --job-name='GMM'
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -o /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/data/21Q1/profiles/GMM/1.fit_gmm.out.txt
#SBATCH -e /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/data/21Q1/profiles/GMM/1.fit_gmm.err.txt
#SBATCH -t 0-06:00:00

cd /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/


ZSCORES="data/21Q1/profiles/CERES_zscores.tsv"
GMM_DIR="data/21Q1/profiles/GMM/"
mkdir -p $GMM_DIR

GMM_OUTPUT="data/21Q1/profiles/GMM/CERES_zscores_2C_with_resp"


source activate superdendrix-env

#python utils/fit_gmm.py -pf ${DIFF_DEP} -o ${GMM_OUTPUT}
python utils/fit_gmm.py -pf ${ZSCORES} -o ${GMM_OUTPUT}
