#!/bin/bash

#SBATCH --partition raphael
#SBATCH --job-name='CT.1'
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -o /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/experiments/20Q2/gmm/logs/pp.1.out.txt
#SBATCH -e /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/experiments/20Q2/gmm/logs/pp.1.err.txt
#SBATCH -t 0-05:00:00

cd /n/fs/ragr-research/projects/sd/gitrepo/superdendrix
#0-764
#0-793
RUNID='gmm'
#RUNID='dependency_probability'
#RUNID='ash'
ALT='MUT'
K=3
CYCLE=10000
RUNDIR='experiments/21Q1/'$RUNID'/cycle'$CYCLE'_BIC0/'
#OUTDIR='experiments/20Q2/'$RUNID'/cycle'$CYCLE'/k'$K'/'$ALT'/'
OUTDIR='experiments/21Q1/'$RUNID'/cycle'$CYCLE'_BIC0/k'$K'/'$ALT'/'
#OUTDIR='experiments/20Q2/'$RUNID'/final/k'$K'/'$ALT'/'
source activate superdendrix-env

RS=1764

#NUMPROFILES=740 # 6sig and BIC 0
NUMPROFILES=765 # 6sig and BIC 0
#NUMPROFILES=511 # 6sig and BIC 10

SWITCH=1
numfiles=$(ls $OUTDIR | wc -l)

if [ $numfiles -eq $NUMPROFILES ];
then
    SWITCH=0
    echo finished
fi

CATFN='allresults_'$ALT'_k'$K'.tsv'
#RUNDIR='experiments/20Q2/'$RUNID'/cycle'$CYCLE'_BIC/'
#RUNDIR='experiments/20Q2/'$RUNID'/cycle'$CYCLE'/'
#RUNDIR='experiments/20Q2/'$RUNID'/final/'

###HEADER="target\taberrations\t#sample\tscores\tz\tmax_score\tz\/max_score(%)\tcoverage\tIC\tranksum_pval\tbrowser_link\tt_opt\tt_total\tP(z>=curve_z)\tP(z<=curve_z)\tP(z>=profile_perm)\tP(z<=profile_perm)\tcurve_avg_rand_z\tcurve_stdev_rand_z\tprofileperm_avg_rand_z\tprofileperm_stdev_rand_zprofile\tfeatures\t#sample\tscores\tW(M)\/max_score(%)\tcoverage\tranksum_pval\tbrowser_link\tt_opt\tt_total\tP_value"

HEADER="profile\tfeatures\t#sample\tscores\tW(M)\/max_score(%)\tcoverage\tranksum_pval\tbrowser_link\tt_opt\tt_total\tP_value\tdirection\tfeature_count"
#HEADER="profile\tfeatures\t#sample\tscores\tW(M)\/max_score(%)\tcoverage\tranksum_pval\tbrowser_link\tt_opt\tt_total\tP_value"

echo $HEADER
if [ $SWITCH -eq 0 ]
then
    echo computing FDR
    cat $OUTDIR* > $RUNDIR$CATFN
    sed -i 's/'$HEADER'//g' $RUNDIR$CATFN
    sed -i '1s/^/'$HEADER'/' $RUNDIR$CATFN
    python utils/compute_FDR.py -i $RUNDIR$CATFN -p 10 -o $RUNDIR'FDR_'$CATFN
else
    echo unfinished
fi

