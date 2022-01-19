#!/bin/bash

#SBATCH --partition raphael
#SBATCH --job-name='CT.1'
#SBATCH -a 2-768
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -o /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/experiments/21Q1/gmm/logs_BIC0/CT.k5.p10000.1.out.%a.txt
#SBATCH -e /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/experiments/21Q1/gmm/logs_BIC0/CT.k5.p10000.1.err.%a.txt
#SBATCH -t 0-12:00:00

cd /n/fs/ragr-research/projects/sd/gitrepo/superdendrix
RUNID='gmm'
ALT='CT'
K=5
CYCLE=10000
CP=10000

#MYID=$((SLURM_ARRAY_TASK_ID * 3))
#END=$((MYID + 2))

EVENTS='data/21Q1/features/IOmerged/CERES_OncoKB_'$ALT'_merged.tsv'	#Ceres alterations
#EVENTS='data/21Q1/features/IOmerged/CERES_full_'$ALT'.tsv'	#Ceres alterations

#MUTLIST='data/21Q1/features/OncoKB_mutation_ct_list.tsv'

OUTCOMES='data/21Q1/profiles/GMM/corrected_CERES_zscores_2C_with_resp_2Cscores.tsv'
#OUTCOMES='data/21Q1/raw/Achilles_gene_dependency.csv'

#OUTDIR='experiments/21Q1/'$RUNID'/cycle'$CYCLE'_full/k'$K'/'$ALT'/'
OUTDIR='experiments/21Q1/'$RUNID'/cycle'$CYCLE'_BIC0/k'$K'/'$ALT'/'
mkdir -p $OUTDIR

NULLMAT='data/21Q1/features/rand_mat/'$ALT'/' #full mut NONSILENT

module load gurobi/9.0.0

#grb_ts
source activate superdendrix-env

GENES='data/21Q1/profiles/TMM/differential_dependencies_tmm.tsv'

RS=1764
G=$(awk -v k=$SLURM_ARRAY_TASK_ID -F $'\t' 'FNR == k {print $1}' $GENES)
D=$(awk -v k=$SLURM_ARRAY_TASK_ID -F $'\t' 'FNR == k {print $2}' $GENES)
#G=$(awk -v k=$SLURM_ARRAY_TASK_ID -F $',' 'FNR == k {print $1}' $GENES)
#D=$(awk -v k=$SLURM_ARRAY_TASK_ID -F $',' 'FNR == k {print $3}' $GENES)
echo $G
echo $D
#python superdendrix/superdendrix_curve_test.py -t 4 -m $EVENTS -T $OUTCOMES -gf $MUTLIST -Tc $G -p $CYCLE -cp $CP -d $D -k $K -nm $NULLMAT -rs $RS -x -curve -o $OUTDIR$G"_cyc"$CYCLE"_k"$K".txt"  # original SD

if [ -f $OUTDIR$G".tsv" ]
then
    echo $OUTDIR$G".tsv EXISTS"
else
    echo $G "not run yet,working now"
    python src/superdendrix.py -t 4 -m $EVENTS -T $OUTCOMES -Tc $G -p $CYCLE -cp $CP -d $D -k $K -nm $NULLMAT -rs $RS -x -curve -o $OUTDIR$G".tsv"  # original SD
    #python src/superdendrix.py -t 4 -m $EVENTS -T $OUTCOMES -Tc $G -p $CYCLE -cp $CP -d $D -k $K -nm $NULLMAT -rs $RS -x -curve -o $OUTDIR$G".tsv"  # original SD
    #python src/superdendrix.py -t 4 -m $EVENTS -T $OUTCOMES -gf $MUTLIST -Tc $G -p $CYCLE -cp $CP -d positive -k $K -nm $NULLMAT -rs $RS -x -curve -o $OUTDIR$G".tsv"  # original SD
fi



#python gitrepo/superdendrix/src/multi-superdendrix.py -t 4 -m $EVENTS -T $OUTCOMES -gf $MUTLIST -Tc $G -p $CYCLE -cp $CP -d $D -k $K -nm $NULLMAT -x -curve -rs $RS -nsub $SUBN -subp $SUBP -o $OUTDIR$G"_cyc"$CYCLE"_k"$K".txt"  # multi superdendrix

#for ((i=$MYID;i<=$END;i++));
#do
#    G=$(awk -v k=$i -F $'\t' 'FNR == k {print $1}' $GENES)
#    echo $G
#    if [ -f $OUTDIR$G"_cyc"$CYCLE"_k"$K".txt" ]
#    then
#        echo $OUTDIR$G"_cyc"$CYCLE"_k"$K".txt EXISTS"
#    else
#        echo $G "not run yet, working now"
#        python superdendrix/superdendrix_curve_test.py -t 4 -m $EVENTS -T $OUTCOMES -Tc $G -p $CYCLE -cp $CP -k $K -nm $NULLMAT -x -curve -o $OUTDIR$G"_cyc"$CYCLE"_k"$K".txt"  # original SD
#        #python old_superdendrix/superdendrix_subopt.py -t 4 -m $EVENTS -T $OUTCOMES -Tc $G -p $CYCLE -cp $CP -k $K -x -o $OUTDIR$G"_cyc"$CYCLE"_k"$K".txt" # all subopt
#    fi
#done
#
###NUMPROFILES=492
####NUMRESULTS=$(ls $OUTDIR* | wc -l)
###
###
###SWITCH=1
###numfiles=$(ls $OUTDIR | wc -l)
###if [ $numfiles -eq $NUMPROFILES ];
###then
###    SWITCH=0
###    echo finished
###fi
###
###CATFN='allresults_'$ALT'.tsv'
###RUNDIR='experiments/superdendrix/'$RUNID'/'
###
###if [ $SWITCH -eq 0 ]
###then
###    echo computing FDR
###    cat /n/fs/ragr-research/projects/sd/experiments/superdendrix/$RUNID/results_$ALT/* > $RUNDIR$CATFN
###    sed -i 's/target\taberrations\t#sample\tscores\tz\tmax_score\tz\/max_score(%)\tcoverage\tIC\tranksum_pval\tbrowser_link\tt_opt\tt_total\tP(z>=curve_z)\tP(z<=curve_z)\tP(z>=profile_perm)\tP(z<=profile_perm)\tcurve_avg_rand_z\tcurve_stdev_rand_z\tprofileperm_avg_rand_z\tprofileperm_stdev_rand_z//g' $RUNDIR$CATFN
###    sed -i '1s/^/target\taberrations\t#sample\tscores\tz\tmax_score\tz\/max_score(%)\tcoverage\tIC\tranksum_pval\tbrowser_link\tt_opt\tt_total\tP(z>=curve_z)\tP(z<=curve_z)\tP(z>=profile_perm)\tP(z<=profile_perm)\tcurve_avg_rand_z\tcurve_stdev_rand_z\tprofileperm_avg_rand_z\tprofileperm_stdev_rand_z/' $RUNDIR$CATFN
###    python superdendrix/compute_FDR.py -i $RUNDIR$CATFN -p 13 -o $RUNDIR'FDR_'$CATFN
###
###else
###    echo unfinished
###fi

#source deactivate
#grb_ts -s
