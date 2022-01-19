#!/bin/bash

#SBATCH --partition raphael
#SBATCH --job-name='mat.MUT.1'
#SBATCH -a 0-999
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -o /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/data/21Q1/features/rand_mat/logs/MUT.OncoKB.merged.out.%a.txt
#SBATCH -e /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/data/21Q1/features/rand_mat/logs/MUT.OncoKB.merged.err.%a.txt
#SBATCH -t 0-01:00:00

cd /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/

EVENTS='data/21Q1/features/IOmerged/CERES_OncoKB_MUT_merged.tsv'	#Ceres MUT FULL CCLE

LOGDIR='data/21Q1/features/rand_mat/logs/'	#Ceres MUT FULL CCLE
mkdir -p ${LOGDIR}
OUTDIR='data/21Q1/features/rand_mat/MUT/'	#Ceres MUT FULL CCLE
mkdir -p ${OUTDIR}

CYCLE=10

PRE=$((SLURM_ARRAY_TASK_ID * 10))
RS=$SLURM_ARRAY_TASK_ID

source activate superdendrix-env
#python superdendrix/generate_null_matrices.py -m $EVENTS -T $OUTCOMES -Tc KRAS -p $CYCLE -curve -o $OUTDIR -pre $PRE -rs $RS
#python superdendrix/generate_null_matrices.py -m $EVENTS -T $OUTCOMES -Tc KRAS -p $CYCLE -o $OUTDIR -pre $PRE -rs $RS
python utils/generate_null_matrices.py -m $EVENTS -p $CYCLE -o $OUTDIR -pre $PRE -rs $RS
