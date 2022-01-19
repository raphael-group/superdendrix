
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# -*- coding: utf-8 -*-

from urllib.parse import quote, urlencode, urlparse
from scipy.stats import ranksums, ttest_ind

import sys, time, argparse, socket, logging
from collections import Counter
from i_o import load_mutation_data, getLogger, load_profiles, load_events
import numpy as np
import math
import random
from collections import defaultdict

# Load our local Python module for computing revealer IC


# Parse arguments
def get_parser():
    description = 'Given mutation profiles and a continous target profiles (and integer k), find the set of genes (of size k) that maximizes weighted target coverage or (-x) the supervised Dendrix score'
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    # data: mutation, profile
    # input: mutation set, profile name, direction

    # Mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True, help='File name for mutation data')
    parser.add_argument('-gf', '--gene_file', default=None, help='File of events to be included (optional).')
    # specific
    parser.add_argument('-q', '--query', required=True, help='query results')
    parser.add_argument('-T', '--target', required=True, help='File name for target profile')
    parser.add_argument('-Tf', '--target_format', default='achilles', choices=['achilles', 'revealer'], help='format of target profile [achilles]')
    parser.add_argument('-Tc', '--target_column', required=False, help='Name of column to use in target profile')
    parser.add_argument('-o', '--output_file', help='Name of output file (csv)')
    # scoring function
    parser.add_argument('-d', '--direction', default='positive', choices=['positive', 'negative'], help='direction: match [positive] or negative values')
    # general
    parser.add_argument('-v', '--verbosity', default=logging.INFO, type=int, help='Flag verbose output')
    parser.add_argument('-rs', '--random_seed', default=1764, type=int, help='Seed for random number generator')
    parser.add_argument('-t', '--threads', type=int, default=-1, help='Number of threads to use ([-1] = all).')

    return parser

args = get_parser().parse_args(sys.argv[1:])
ind_profile = 1
ind_feature = 2
ind_dir = 12

# load data
with open(args.query) as f:
    head=f.readline()
    arrs = [l.rstrip().split("\t") for l in f if not l.startswith("target")]


mutationMatrix = args.mutation_matrix

random.seed(args.random_seed) # founding year of Brown University
np.random.seed(args.random_seed)

# Generate/load the target profile
eventToCases, mutation_samples = load_events(args.mutation_matrix, verbose=1)
profiles, profile_samples = load_profiles(args.target, sample_whitelist=mutation_samples,verbose=1)

# Load the mutation data
mutations = load_mutation_data(mutationMatrix, profile_samples, args.gene_file)#, args.mutations_only, args.min_freq, args.max_freq)
m, N, features, samples, events_to_samples, samples_to_genes = mutations


#newprofiles = profiles

newprofiles = dict()

for k in profiles.keys():
    newk = k.split("_")[0]
    newprofiles[newk] = profiles[k]

# start checking

def run(target, module, direction, newprofiles):
    if "_" in target:
        gene = target.split("_")[0]
    else:
        gene = target
    result_curr = [target,str(module),direction]

    coef = 1
    if direction == 'negative': coef = -1

    if gene not in newprofiles.keys():
        print("profile is not in the data")
        result_curr.append("FALSE")
    
    if gene in newprofiles.keys():
        profile = dict((s, coef * newprofiles[gene][i]) for i, s in enumerate(profile_samples) if s in samples)
        result_curr.append("TRUE")

    print("-------------------")
    print("curr gene: ",gene)
    print("curr feature set: ",module)
    print("curr direction: ",direction)
    print('\n')
    print("check each feature is in the feature data")
    isIn = []
    for i in module:
        #print("module %s is in features? " % i,i in features)
        isIn.append(i+":"+str(i in features))
    result_curr.append(str(isIn))

    print("\n")

    if gene not in newprofiles.keys():
        result_curr.append("NA")
        return(result_curr)

    # we multiply by minus one because the weights had been transformed that way (effect is just changing the sign of IC, but this way it is correct)



    ##Wilcoxon rank sum test and two-sample t-test
    if len(module) == 0:
        ranksum_pval = 1
    else:
        mutsamples = set()
        for i in range(len(module)):
            mutsamples = mutsamples.union(events_to_samples[module[i]])
        nomutsamples = set(samples) - mutsamples

        mutscores = [profile[s] for s in mutsamples]
        nomutscores = [profile[s] for s in nomutsamples]

        ranksum_pval = ranksums(mutscores, nomutscores)[1]
        t_pval = ttest_ind(mutscores, nomutscores, equal_var = False, nan_policy='omit')[1]


    result_curr.append(str(ranksum_pval))
    result_curr.append(str(t_pval))
    return(result_curr)

results = []
for arr in arrs:
    target = arr[ind_profile]
    #gene = target.split("_")[0]
    module = arr[ind_feature].split(",")
    direction = arr[ind_dir]
    result_line = run(target,module,direction,newprofiles)
    results.append(result_line)



##output
if args.output_file:
    of = open(args.output_file, 'w')
    of.write("profile\tfeatures\tdirection\tprofile_is_in_data\tfeature_is_in_data\tranksum_pval\tt_pval\n")
    for curr in results:
        of.write("\t".join(curr) + "\n")
