#!/usr/bin/env python

import sys, time, argparse, os
import statistics

#sys.path.append(os.path.join(os.path.dirname(__file__), '../../utils'))
from i_o import *

from scipy.stats import ranksums, ttest_ind

# parse SuperDendrix results 
# input: tsv file from SuperDendrix
# output: 
# 1) list of significant profiles,
# 2) list of significant mutations
# 3) pairs of associations from significant associations

# Parse arguments
def get_parser():
    description = 'Parse SuperDendrix results'
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    parser.add_argument('-p', '--profiles', required=True, help='File name for phenotype')
    # Mutation data
    parser.add_argument('-i', '--results', required=True, help='File name for SD results')
    parser.add_argument('-ef', '--mutation_matrix', required=False, help='File name for mutation data')
    parser.add_argument('-o', '--output_file', required=False, help='Directory name for output')
    parser.add_argument('-ct', '--CT', default=False,action='store_true')
    parser.add_argument('-find', '--fdr_index', default=13, required=False, action='store_true')
    return parser

args = get_parser().parse_args(sys.argv[1:])

def load_single_profile(args, gene):
    # Generate/load the target profile
    w = {} # weights
    global numNanInf
    numNanInf = 0
    with open(args.profiles) as f:
        line = f.readline()
        delim = "\t" if args.profiles.endswith(".tsv") else ","
        ind = line.rstrip().split(delim).index(gene)# + 1
        arrs = [ l.rstrip().split(delim) for l in f if not l.startswith("#") ]
        #arrs = [ re.findall(r"[-\w']+", l) for l in f if not l.startswith("#") ]
        #for arr in [arr for arr in arrs if arr[0] in patients]:
        for arr in arrs:
            if arr[ind] == "":
                w[arr[0]] = 0.0
            else:
                w[arr[0]] = float(arr[ind])
    return w

if args.mutation_matrix:
    m, N, features, samples, features_to_samples, samples_to_features = load_mutation_data(args.mutation_matrix)

pairs = []
profile_genes = []
mutated_genes = dict()
mutated_features = dict()

rows = []
with open(args.results, 'r') as f:
    if args.results.endswith("tsv"):
        delim = "\t"
    else:
        delim = ","
    arrs = [l.rstrip('\n').split(delim) for l in f]
    header = arrs[0]
    rest = arrs[1:]
    for arr in rest:
        target = arr[0]
        print(target)
        profile = load_single_profile(args, target)
        feature1 = arr[1]
        feature2 = arr[2]
        cl1 = features_to_samples[feature1]
        cl2 = features_to_samples[feature2]
        cl_union = cl1.union(cl2)

        m1_score = sum([ profile[sample] for sample in cl1 ])
        m2_score = sum([ profile[sample] for sample in cl2 ])
        total_score = sum([ profile[sample] for sample in cl_union ])
        m1_contribution = m1_score/total_score
        m2_contribution = m2_score/total_score
        arr = arr + [str(m1_contribution), str(m2_contribution)]
        rows.append(arr)

# write results
# let's write mutations and cancer types separately
if args.output_file:
    print("writing files now")
    of = open(args.output_file + "pairs.tsv", 'w')
    of.write('\t'.join(header) + "\n")
    for row in rows:
        of.write('\t'.join(row) + "\n")
    of.close()
