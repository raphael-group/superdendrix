#!/usr/bin/env bash
cd /n/fs/ragr-research/projects/sd/gitrepo/oncokb-annotator

TOKEN="12cd9caf-9fb8-4098-80fd-9de027fe6b9a" #OncoKB API Token

IMAF="../data/raw/CCLE_mutations_avana.tsv"
OMAF="../data/raw/CCLE_mutations_avana_oncoKB.tsv"

echo "starting"
python ../oncokb-annotator/MafAnnotator.py -i $IMAF -o $OMAF -b $TOKEN -q HGVSp_Short
