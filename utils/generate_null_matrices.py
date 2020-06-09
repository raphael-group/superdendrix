#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, time, argparse, socket, logging
from collections import Counter
from i_o import load_mutation_data, getLogger, load_profiles, load_events
import numpy as np
import math
import random
from collections import defaultdict
from curveball_main import *
from curveball import *
from scanutils import *

# Parse arguments
def get_parser():
    description = 'Given mutation profiles and a continous target profiles (and integer k), find the set of genes (of size k) that maximizes weighted target coverage or (-x) the supervised Dendrix score'
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    # Mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True, help='File name for mutation data')
    parser.add_argument('-min_f', '--min_freq', type=int, default=0, help='Minimum gene mutation frequency [0].')
    parser.add_argument('-max_f', '--max_freq', type=int, default=float("inf"), help='Maximum gene mutation frequency [infty].')
    parser.add_argument('-pf', '--sample_file', default=None, help='File of samples to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None, help='File of events to be included (optional).')
    # specific
    parser.add_argument('-M', '--mutations_only', default=False, action="store_true", help='Consider only mutations (events ending with \'_MUT\')')
    parser.add_argument('-o', '--output_file', help='Name of output file (csv)')
    parser.add_argument('-p', '--p_val_it', type=int, default=-1, help='Number of iterations to compute empirical p-value ([-1] = don\'t)')
    parser.add_argument('-pre', '--prefix', type=int, default=0, help='prefix of output file')
    # general
    parser.add_argument('-v', '--verbosity', default=logging.INFO, type=int, help='Flag verbose output')
    parser.add_argument('-rs', '--random_seed', default=1764, type=int, help='Seed for random number generator')
    return parser

args = get_parser().parse_args(sys.argv[1:])
t_init = time.time()
# Set up logger
logger = getLogger(args.verbosity)
logger.info("# calling %s" % " ".join(sys.argv))
logger.info("# at time %s" % time.strftime("%H:%M:%S on %a, %d %b %Y "))
logger.info("# on machine %s" % socket.gethostname())

# 
global k, N, events_to_samples, samples_to_genes, samples, module
mutationMatrix = args.mutation_matrix
seed = []

random.seed(args.random_seed) # founding year of Brown University
np.random.seed(args.random_seed)


# Load the mutation data
m, N, genes, samples, events_to_samples, samples_to_genes = load_mutation_data(mutationMatrix, None, args.gene_file, args.mutations_only, args.min_freq, args.max_freq)

logger.info('* Mutation data')
logger.info('\t- Alterations: %s' % m)
logger.info('\t- Samples: %s' % N)
logger.info('* Seed genes: %s' %  ','.join(seed))

# matrix permutation using curveball method
t_s = time.time()
mat = np.empty((len(samples),len(genes)))
print(len(samples),len(genes))

# convert to matrix
for i in range(len(samples)):
    row = [mut in samples_to_genes[samples[i]] for mut in genes]
    mat[i] = row

# input numpy matrix, find matrices
presence = find_presences(mat)
t_ss = time.time()
for i in range(args.p_val_it):
    permuted_mat = curve_ball(mat, presence)
    permuted_samples_to_genes = []
    for r in range(len(permuted_mat)):
        sample = samples[r]
        mutations = [s for i, s in zip(permuted_mat[r], genes) if i == 1]
        row = [sample] + mutations
        permuted_samples_to_genes.append(row)

    of = open(args.output_file + str(i + args.prefix) + ".tsv", 'w')
    of.write("#Samples\tEvents\n")
    for row in permuted_samples_to_genes:
        of.write("\t".join(row) + "\n")
    of.close()
    print(i)
    if i==100:
        print(str((time.time() - t_ss)/60),"min for 100 cycles")

print(str((time.time() - t_s)/60),"min")
