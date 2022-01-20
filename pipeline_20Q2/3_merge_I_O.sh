#!/usr/bin/env bash

cd /n/fs/ragr-research/projects/sd/gitrepo/superdendrix/
'''
FEATUREDIR="data/21Q1/features/"
ODIR=$FEATUREDIR"IOmerged/"
mkdir -p $ODIR

IMAF="data/21Q1/features/CERES_OncoKB_MUT.tsv"
OMAF="data/21Q1/features/IOmerged/CERES_OncoKB_MUT_merged.tsv"

echo "starting"
python utils/merge_I_O.py $IMAF $OMAF

IMAF="data/21Q1/features/CERES_OncoKB_CT.tsv"
OMAF="data/21Q1/features/IOmerged/CERES_OncoKB_CT_merged.tsv"

echo "starting"
python utils/merge_I_O.py $IMAF $OMAF


IMAF="data/21Q1/features/CERES_full_MUT.tsv"
OMAF="data/21Q1/features/IOmerged/CERES_full_MUT_merged.tsv"

echo "starting"
python utils/merge_I_O.py $IMAF $OMAF


IMAF="data/21Q1/features/CERES_full_CT.tsv"
OMAF="data/21Q1/features/IOmerged/CERES_full_CT_merged.tsv"

echo "starting"
python utils/merge_I_O.py $IMAF $OMAF
'''

IMAF="data/21Q1/features/CERES_full_CT_lung.tsv"
OMAF="data/21Q1/features/IOmerged/CERES_full_CT_lung_merged.tsv"

echo "starting"
python utils/merge_I_O.py $IMAF $OMAF

IMAF="data/21Q1/features/CERES_OncoKB_CT_lung.tsv"
OMAF="data/21Q1/features/IOmerged/CERES_OncoKB_CT_lung_merged.tsv"

echo "starting"
python utils/merge_I_O.py $IMAF $OMAF


#python FusionAnnotator.py -i $IF -o $OF -c $IC -b $TOKEN
#python CnaAnnotator.py -i $ICNA -o $OCNA -c $IC -b $TOKEN
#python ClinicalDataAnnotator.py -i $IC -o $OC -a $OMAF,$OCNA,$OF
#python OncoKBPlots.py -i $OC -o $OCPDF -c ONCOTREE_CODE #-n 10
#python GenerateReadMe.py -o $README
