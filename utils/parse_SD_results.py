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
    # Mutation data
    parser.add_argument('-i', '--results', required=True, help='File name for SD results')
    parser.add_argument('-ef', '--mutation_matrix', required=False, help='File name for mutation data')
    parser.add_argument('-o', '--output_file', required=False, help='Directory name for output')
    parser.add_argument('-ct', '--CT', default=False,action='store_true')
    parser.add_argument('-find', '--fdr_index', default=13, required=False, action='store_true')
    return parser

args = get_parser().parse_args(sys.argv[1:])

if args.mutation_matrix:
    m, N, features, samples, features_to_samples, samples_to_features = load_mutation_data(args.mutation_matrix)


pairs = []
profile_genes = []
mutated_genes = dict()
mutated_features = dict()
with open(args.results, 'r') as f:
    if args.results.endswith("tsv"):
        delim = "\t"
    else:
        delim = ","
    arrs = [l.rstrip('\n').split(delim) for l in f if not l.startswith('profile')]
    for arr in arrs:
        profile = arr[0]
        features = arr[1].split(",")
        #print(arr[args.fdr_index])
        fdr = float(arr[args.fdr_index])
        if fdr <= 0.20:
            profile_genes.append([profile,str(fdr)])
            for f in features:
                if args.CT:
                    if f.endswith("_Cancer"):
                        mut_gene, mut_class = f.split("_Cancer")
                else:
                    if f.endswith("_MUT"):
                        mut_gene, mut_class = f.split("_MUT")
                if mut_gene in mutated_genes:
                 mutated_genes[mut_gene] += 1
                else:
                    mutated_genes[mut_gene] = 1
                if f in mutated_features:
                    mutated_features[f] += 1
                else:
                    mutated_features[f] = 1
                pairs.append([profile,f, str(len(features_to_samples[f])), str(fdr)])

#            if args.CT:
#                if f.endswith("Cancer"):
#                    pairs.append([profile,f,profile+","+f,profile+","+f])
#            else:
#                pairs.append([profile,f,profile+","+f,profile+","+f[:-6], str(len(features_to_samples[f]))])



# write results
# let's write mutations and cancer types separately
if args.output_file:
    print("writing files now")
    of = open(args.output_file + "pairs.tsv", 'w')
    header = 'profile\tfeature\tcoverage\tFDR\n'
    of.write(header)
    for row in pairs:
        of.write('\t'.join(row) + "\n")
    of.close()
    
    of1 = open(args.output_file + "profiles.tsv", "w")
    header = "profile\tfdr\n"
    of1.write(header)
    for row in profile_genes:
        of1.write("\t".join(row) + "\n")
    of1.close()

    of2 = open(args.output_file + "mutated_genes.tsv", "w")
    header = "mutated_gene\tcount\n"
    of2.write(header)
    for k,v in mutated_genes.items():
        of2.write(k + "\t" + str(v) + "\n")
    of2.close()

    of3 = open(args.output_file + "mutated_features.tsv", "w")
    header = "feature\tcount\ttype\tgene\n"
    of3.write(header)
    for k,v in mutated_features.items():
        if k.endswith("_Cancer"):
            gene, feature_type = k.split("_Cancer")
            of3.write(k + "\t" + str(v) + "\t" + feature_type + "\t" + gene + "\n")
    of3.close()


