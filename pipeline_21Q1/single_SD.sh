#!/bin/bash

cd /n/fs/ragr-research/projects/sd/gitrepo/superdendrix
ALT='MUT'
K=3
CYCLE=10000
CP=10000

EVENTS='data/21Q1/features/IOmerged/CERES_OncoKB_'$ALT'_merged.tsv'	#Ceres alterations

OUTCOMES='data/21Q1/profiles/GMM/corrected_CERES_zscores_2C_with_resp_2Cscores.tsv'


OUTDIR='experiments/21Q1/gmm/scrap/'
mkdir -p $OUTDIR

NULLMAT='data/21Q1/features/rand_mat/'$ALT'/' #full mut NONSILENT

#grb_ts
source activate superdendrix-env

GENES='data/21Q1/profiles/TMM/differential_dependencies_tmm.tsv'

RS=1764
#G="MYB_(4602)"
G='USP18_(11274)'
D="negative"
python src/superdendrix.py -t 4 -m $EVENTS -T $OUTCOMES -Tc $G -p $CYCLE -cp $CP -d $D -k $K -nm $NULLMAT -rs $RS -x -curve -o $OUTDIR$G".tsv"

#for i in {2..512..1}
#do
#    G=$(awk -v k=$i -F $'\t' 'FNR == k {print $1}' $GENES)
#    D=$(awk -v k=$i -F $'\t' 'FNR == k {print $2}' $GENES)
#    echo $G
#    python src/superdendrix.py -t 4 -m $EVENTS -T $OUTCOMES -gf $MUTLIST -Tc $G -p $CYCLE -cp $CP -d $D -k $K -nm $NULLMAT -rs $RS -x -curve -o $OUTDIR$G".tsv"
#done

