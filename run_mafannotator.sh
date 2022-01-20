#!/usr/bin/env bash

#cd /n/fs/ragr-research/projects/sd/gitrepo/oncokb-annotator

TOKEN="12cd9caf-9fb8-4098-80fd-9de027fe6b9a" #OncoKB API Token

#IMAF="../superdendrix/data/20Q4/raw/CCLE_mutations.csv"
#OMAF="../superdendrix/data/20Q4/features/CCLE_mutations_20Q4_oncoKB.tsv"
#IMAF="../superdendrix/data/integrated/raw/CCLE_mutations_merged.csv"
#OMAF="../superdendrix/data/integrated/features/CCLE_mutations_integrated_oncoKB.tsv"

IMAF="data/raw/CCLE_mutations_avana.tsv"
#OMAF="../superdendrix/data/20Q2/features/CCLE_mutations_oncoKB_3.tsv"
OMAF="data/raw/CCLE_mutations_avana_oncoKB.tsv"



##IMAF="../superdendrix/data/20Q4/features/NFE2L2_mutations_20Q4.csv"
##OMAF="../superdendrix/data/20Q4/features/NFE2L2_mutations_20Q4_oncoKB_38.tsv"
README="data/example_README.txt"

mkdir -p OMAF

echo "starting"
#python ../../oncokb-annotator/MafAnnotator.py -i $IMAF -o $OMAF -b $TOKEN -q HGVSp_Short
python oncokb-annotator/MafAnnotator.py -i $IMAF -o $OMAF -b $TOKEN -q HGVSp_Short

#python FusionAnnotator.py -i $IF -o $OF -c $IC -b $TOKEN
#python CnaAnnotator.py -i $ICNA -o $OCNA -c $IC -b $TOKEN
#python ClinicalDataAnnotator.py -i $IC -o $OC -a $OMAF,$OCNA,$OF
#python OncoKBPlots.py -i $OC -o $OCPDF -c ONCOTREE_CODE #-n 10
#python GenerateReadMe.py -o $README
