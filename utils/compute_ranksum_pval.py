
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# -*- coding: utf-8 -*-

from urllib.parse import quote, urlencode, urlparse
from scipy.stats import ranksums

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
    parser.add_argument('-q', '--query', required=True, help='query alteration')
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

gene = args.target_column.split("_")[0]
print("------------")
print("profile is: ",gene)




t_init = time.time()
# Set up logger
  
logger = getLogger(args.verbosity)
#logger.info("# calling %s" % " ".join(sys.argv))
#logger.info("# at time %s" % time.strftime("%H:%M:%S on %a, %d %b %Y "))
#logger.info("# on machine %s" % socket.gethostname())
    
    # 
global N, events_to_samples, samples_to_genes, samples, module
mutationMatrix = args.mutation_matrix
    
random.seed(args.random_seed) # founding year of Brown University
np.random.seed(args.random_seed)

# Generate/load the target profile

eventToCases, mutation_samples = load_events(args.mutation_matrix, verbose=1)
profiles, profile_samples = load_profiles(args.target, sample_whitelist=mutation_samples,verbose=1)

# Load the mutation data

mutations = load_mutation_data(mutationMatrix, profile_samples, args.gene_file)#, args.mutations_only, args.min_freq, args.max_freq)
m, N, features, samples, events_to_samples, samples_to_genes = mutations
coef = 1
if args.direction == 'negative': coef = -1

global profile

newprofiles = dict()
for k in profiles.keys():
    newk = k.split("_")[0]
    newprofiles[newk] = profiles[k]

profile = dict((s, coef * newprofiles[args.target_column][i]) for i, s in enumerate(profile_samples) if s in samples)
    
logger.info('* Mutation data')
logger.info('\t- Alterations: %s' % m)
logger.info('\t- Samples: %s' % N)
    

t_start = time.time()

module = args.query.split(",")

print("\n")
print("check each feature is in the feature data")
isIn = []
for i in module:
    print("module %s is in features? " % i,i in features)
    isIn.append(i+":"+str(i in features))
cases_by_event = [events_to_samples[event] for event in module]
print("feature set is :",module)
print("\n")

mut_samples = { s for g in module for s in events_to_samples[g] }
# we multiply by minus one because the weights had been transformed that way (effect is just changing the sign of IC, but this way it is correct)

##Wilcoxon rank sum test
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

logger.info('Ranksum_pval: %f' % ranksum_pval)
##output
if args.output_file:
    of = open(args.output_file, 'w')
    of.write("profile\tfeatures\tranksum_pval\tdirection\tfeature_is_in_data\n")
    of.write(str(args.target_column) + "\t")
    of.write(",".join(module))
    of.write("\t" + str(ranksum_pval))
    of.write("\t" + str(args.direction))
    of.write(",".join(isIn))

    
   
for g in module:
    logger.info('%s %s' % (g, len(events_to_samples[g])))
