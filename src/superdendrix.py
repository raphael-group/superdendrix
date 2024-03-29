#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Supervised Dendrix model for explaining the positive parts of a weighted target profile with mutually exclusive events.
GW Klau, M Leiserson, B Raphael
"""
__author__ = """Gunnar Klau (gunnar.klau@cwi.nl)"""
__date__ = ""
__credits__ = ""
__revision__ = ""
#    Copyright (C) 2016 by
#    Gunnar Klau (gunnar.klau@cwi.nl)
#    All rights reserved.
#    MIT license.

from scipy.stats import ranksums

import sys, time, argparse, socket, logging, os
from collections import Counter
from gurobipy import *
import numpy as np
import math
import random
import statistics
from collections import defaultdict

from curveball_main import *
from curveball import *
from scanutils import *

from i_o import load_mutation_data, getLogger, load_profiles, load_events

# Parse arguments
def get_parser():
    description = 'Given mutation profiles and a continous target profiles (and integer k), find the set of genes (of size k) that maximizes weighted target coverage or (-x) the supervised Dendrix score'
    parser = argparse.ArgumentParser(description=description, fromfile_prefix_chars='@')
    # Mutation data
    parser.add_argument('-m', '--mutation_matrix', required=True, help='File name for mutation data')
    parser.add_argument('-nm', '--null_matrices', required=False, help='directory for curveball matrices')
    parser.add_argument('-min_f', '--min_freq', type=int, default=0, help='Minimum gene mutation frequency [0].')
    parser.add_argument('-max_f', '--max_freq', type=int, default=float("inf"), help='Maximum gene mutation frequency [infty].')
    parser.add_argument('-pf', '--sample_file', default=None, help='File of samples to be included (optional).')
    parser.add_argument('-gf', '--gene_file', default=None, help='File of events to be included (optional).')
    # specific
    parser.add_argument('-k', '--cardinality', default=-1, type=int, help="Fixed cardinality of event subset [-1]")
    parser.add_argument('-M', '--mutations_only', default=False, action="store_true", help='Consider only mutations (events ending with \'_MUT\')')
    parser.add_argument('-T', '--target', required=True, help='File name for target profile')
    parser.add_argument('-Tc', '--target_column', required=False,type=str, help='Name of column to use in target profile')
    parser.add_argument('-o', '--output_file', help='Name of output file (csv)')
    parser.add_argument('-r', '--reoptimize', default=False, action="store_true", help='Reoptimize for number of genes if k=-1')
    parser.add_argument('-cp', '--cond_it', type=int, default=10000, help='Number of iterations to compute conditional p-value ([-1] = don\'t)')
    parser.add_argument('-p', '--p_val_it', type=int, default=10000, help='Number of iterations to compute empirical p-value ([-1] = don\'t)')
    # scoring function
    parser.add_argument('-d', '--direction', default='increased_dependency', choices=['increased_dependency','decreased_dependency','positive', 'negative'], help='direction: match [positive] or negative values')
    parser.add_argument('-u', '--unit_weights', default=False, action="store_true", help='Use unit weights and simple binarization [False]')
    parser.add_argument('-O', '--offset', type=float, default=0, help='manual weight offset [0]')
    # general
    parser.add_argument('-v', '--verbosity', default=logging.INFO, type=int, help='Flag verbose output')
    parser.add_argument('-rs', '--random_seed', default=1764, type=int, help='Seed for random number generator') # founding year of Brown University
    parser.add_argument('-t', '--threads', type=int, default=-1, help='Number of threads to use ([-1] = all).')

    return parser

def run(args):
    t_init = time.time()
    # Set up logger
    logger = getLogger(args.verbosity)
    logger.info("# calling %s" % " ".join(sys.argv))
    logger.info("# at time %s" % time.strftime("%H:%M:%S on %a, %d %b %Y "))
    logger.info("# on machine %s" % socket.gethostname())
    # 
    global k, N, samples, module
    k = args.cardinality
    mutationMatrix = args.mutation_matrix
    
    random.seed(args.random_seed)
    np.random.seed(args.random_seed)
    
    # Generate/load the target profile
    w = {} # weights
    global numNanInf
    numNanInf = 0

    assert(args.target_column)
    with open(args.target) as f:
        line = f.readline()
        delim = "\t" if args.target.endswith(".tsv") else ","
        ind = line.rstrip().split(delim).index(args.target_column)# + 1
        arrs = [ l.rstrip().split(delim) for l in f if not l.startswith("#") ]
        for arr in arrs:
            try:
                if math.isnan(float(arr[ind])):
                    logger.info("Warning: profile of sample %s is NaN: removing sample." % arr[0])
                    numNanInf += 1
                elif math.isinf(float(arr[ind])):
                    w[arr[0]] = 100.0 if float(arr[ind])==float("inf") else -100.0
                    if (args.direction == 'negative' or args.direction == 'increased_dependency'):
                        w[arr[0]] = -w[arr[0]]
                    numNanInf += 1
                elif arr[ind] == "":
                    logger.info("Warning: profile of sample %s is blank: removing sample." % arr[0])
                    numNanInf += 1
                else:
                    w[arr[0]] = float(arr[ind]) + args.offset
                    if (args.direction == 'negative' or args.direction == 'increased_dependency'):
                        w[arr[0]] = -w[arr[0]]
                    if args.unit_weights:
                        if w[arr[0]] <= 0: w[arr[0]] = -1
                        else: w[arr[0]] = +1
            except:
                logger.warning("Warning: %s %s" % (arr[0], arr[ind]))
                numNanInf += 1
    
    # terminate if too many inf or NaN
    if float(numNanInf)/len(arrs) > 0.1:
        print("too many Nan or Inf, terminating")
        logger.info('Nan, Inf fraction: %s ' % (float(numNanInf)/len(arrs)))
        if args.output_file:
            of = open(args.output_file, 'w')
            of.write(str(args.target_column) + "\t")
            of.write("incomplete due to high Nan, Inf\t")
            of.write("Nan, Inf fraction: %s" % (float(numNanInf)/len(arrs)))
        return 1

    profile_samples = w.keys()
    eventToCases, mutation_samples = load_events(args.mutation_matrix, verbose=1)   
    
    # Load the mutation data
    mutations = load_mutation_data(mutationMatrix, profile_samples, args.gene_file, args.mutations_only, args.min_freq, args.max_freq)
    m, N, genes, samples, events_to_samples, samples_to_genes = mutations

    coef = 1
    if (args.direction == 'negative' or args.direction == 'increased_dependency'):
        coef = -1

    profile = dict((s, (w[s])) for i, s in enumerate(samples) if s in samples)
    
    logger.info('* Mutation data')
    logger.info('\t- Alterations: %s' % m)
    logger.info('\t- Samples: %s' % N)
   
    t_start = time.time()

    z, module, best_single, cov, act_cov, act_sample,x,y,mo = optimize_model(genes, samples, profile, events_to_samples, samples_to_genes, args, logger)
    t_opt = str(time.time() - t_start)

    # pruning module
    t_pruning = time.time()
    logger.info("Before model selection, module: " + str(module))

    if args.cond_it != -1:
        while len(module) > 1:
            logger.info("Number of features in module:" + str(len(module)))
            worst_event, worst_count = conditional_permutation_test(module, events_to_samples, profile, samples, args.cond_it)
            #if args.verbosity: print('worst',worst_event, worst_count/float(args.cond_it))
            logger.info('Least significant feature: ' + str(worst_event) + str(worst_count/float(args.cond_it)))
            if worst_count/float(args.cond_it) <= (1.0/float(args.cond_it)): break
            module.remove(worst_event)
    
    logger.info('time for model selection(sec): ' + str((time.time() - t_pruning)))
    logger.info('module after module selection: ' + str(module))

    cases_by_event = [events_to_samples[event] for event in module]
    zP = SuperW(cases_by_event, profile, samples)

    # init for permutation tests
    p_z_arr_curve = []
    p_z_arr_profile = []
    p_val_str_curve = "na"
    #p_val_str_curve_2 = "na"
    
    # matrix permutation using curveball method
    k = len(module)
    t_s1 = time.time()
    if len(module) > 0 and args.p_val_it > 0:
        nullmat_dir = args.null_matrices
        c = 0
        t_curve_arr = []
        t_opt_arr = []
        no_better = 0
        no_worse = 0
   
        for i in range(args.p_val_it):
            t_s = time.time()
            
            p_mutations = load_mutation_data(nullmat_dir+str(i)+".tsv", profile_samples, args.gene_file, args.mutations_only, args.min_freq, args.max_freq)
            p_m, p_N, p_genes, p_samples, p_events_to_samples, p_samples_to_genes = p_mutations
            t_curve = time.time() - t_s

            t_s = time.time()
            p_z, p_module, p_best_single, p_cov, p_act_cov, p_act_sample,p_x,p_y,p_mo = optimize_model(p_genes, p_samples, w, p_events_to_samples, p_samples_to_genes, args, logger) #change args.k
            if p_z >= zP: no_better = no_better + 1
            if p_z <= zP: no_worse = no_worse + 1
            p_z_arr_curve.append(p_z)
        p_val_str_curve = str(float(no_better)/args.p_val_it)
        #p_val_str_curve_2 = str(float(no_worse)/args.p_val_it)
        print("Time for permutation test:",str((time.time() - t_s1)/60),"min")

    logger.info('coverage: %d/%d' % (cov, N))
    if cov > 0: obj_norm = z/cov
    else: obj_norm = 0
    module.sort(key=lambda g: len(events_to_samples[g]), reverse=True)
    max_score = sum(profile[p] for p in samples if profile[p] > 0)
    if max_score == 0:
        max_score = -0.1
    
    t_model = time.time() - t_start
    
    mut_samples = { s for g in module for s in events_to_samples[g] }
    ordered_w   = [ -profile[s] for s in samples ]
    

    ##for individual aberration scores
    global scores_ind
    scores_ind = []
    scores_dict = dict()
    print("-------------------------")
    print(module)
    for i in range(len(module)):
        sampleset = [events_to_samples[module[i]]]
        scores_ind.append(float(round(SuperW(sampleset,profile,samples),2)))
        print("here")
        print(float(round(SuperW(sampleset,profile,samples),2)))
        scores_dict[module[i]] = float(scores_ind[i])
    scores_ind.sort(reverse=True)
    scores_ind = map(str,scores_ind)
    module.sort(key=lambda g: float(scores_dict[g]), reverse=True)
    
    ##for coverage
    cover = set()
    cov_ind = []
    for i in range(len(module)):
        cov_ind.append(str(len(events_to_samples[module[i]])))
    
    for i in range(len(module)):
        cover = cover.union(events_to_samples[module[i]])
    
    #logger.info(args.target_column + "\t" + str(z) + "\t" + str(max_score)  + "\t" + str(z/(max_score+0.0001)) + "\t" + str(cov) + "\t" + str(len(samples)) + "\t" + str(act_cov) + "\t" + str(act_sample) + "\t" + str(obj_norm) + "\t" + str(t_model) + "\t" + str(module) + "\t" + str(len(module)) + "\t" + str(best_single))
    
    ##output
    if args.output_file:
        of = open(args.output_file, 'w')
        of.write("profile\tfeatures\t#samples\tscores\tW(M)(%)\tcoverage\tP_value\tdirection\tfeature_count\n")
        of.write(str(args.target_column) + "\t")
        if len(module) > 0:
            of.write(",".join(module))
            of.write("\t" + ",".join(cov_ind))
            of.write("\t" + ",".join(scores_ind))
        else:
            of.write("nan")
            of.write("\t" + "nan")
            of.write('\t' + "nan")
        of.write("\t" + str(round(100*zP/max_score,4)))
        of.write("\t" + str(len(cover))+"/"+str(N))
        of.write("\t" + str(p_val_str_curve))
        of.write("\t" + str(args.direction))
        of.write("\t" + str(len(module))+"\n")
    
    logger.info("Profile: " + args.target_column)
    logger.info("FINAL MODULE: " + ", ".join(module))
    for g in module:
        logger.info('%s %s' % (g, len(events_to_samples[g])))

    return module


def optimize_model(genes, samples, w, events_to_samples, samples_to_genes, args, logger=getLogger()):
    """Builds and optimizes the model using Gurobi.
      * **bound** (*bool*) - use a bound (opt value has to be better than that bound)
      * **z_bound** (*float*) - bound to use in that case

    **Returns:**
      * **z** (*float*) - optimal solution value.
      * **module** (*list*) - an optimal solution.
      * **cov** (*int*) - no of samples covered by solution
      * **act_cov** (*int*) - no of active samples covered by solution
      * **act_sample** (*int*) - no of total active samples
    """
    try:
#        logger.info('* Finding best event...')
        z_single = float("-inf")
        for e in genes:
            z_e = sum(w[s] for s in events_to_samples[e])
            if z_e > z_single:
                z_single = z_e
                e_single = e

#        logger.info('* Building ILP model...')
        m = Model('superdendrix')
        m.setParam('OutputFlag', int(logger.getEffectiveLevel() > logging.INFO))
        if args.threads != -1: m.params.threads = args.threads
        x, y = {}, {}
        # variables
        for g in genes: x[g] = m.addVar(vtype=GRB.BINARY, name = 'x_%s' % g)
        for p in samples: y[p] = m.addVar(vtype=GRB.BINARY, name = 'y_%s' % p)
        m.update()

        #obj = quicksum(y[p] for p in one_samples) - quicksum(y[p] for p in zero_samples)
        #if k != -1:
        m.setObjective(build_obj(x, y, w, samples_to_genes, args), GRB.MAXIMIZE)
        #else:
        #    m.setObjective(obj - quicksum(x[g] for g in genes)/len(genes), GRB.MAXIMIZE)
        # set the bound constraint if appropriate

        # cardinality constraint
        if k != -1: m.addConstr(quicksum(x[g] for g in genes) <= k)
        # m.addConstr(quicksum(x[g] for g in genes) >= 1)  # option1 for 0 len modules

        # coverage constraints
        for p in samples:
            m.addConstr(y[p] <= quicksum(x[g] for g in samples_to_genes[p]))
            for g in samples_to_genes[p]:
                if w[p] < 0: m.addConstr(y[p] >= x[g])

        m.update()
        #m.write('opt_REVEALER.lp')

        # if verbose: print('* Tuning model...')
        # m.tune()
        # m.getTuneResult(0);
        # m.write('tune.prm')

#        logger.info('* Optimizing model...')
        m.optimize()

        if m.SolCount == 0:
            logger.warning('%s No solution found, optimization status = %d' % (args.target_column, m.Status))
            return float("-inf"), [], "", -1, -1, -1, "n/a", "n/a" # z, module, e_single, cov, act_cov, act_sample, p_val_str, p_val_single_str
        else:
            z = m.ObjVal

            if k == -1 and args.reoptimize: #second optimization phase: min # genes
                logger.info('* Re-optimizing model...')
                m.addConstr(obj == m.ObjVal)
                m.setObjective(quicksum(x[g] for g in genes), GRB.MINIMIZE)
                m.optimize()


        #print m.RunTime, t_ILP
        ## output module ##
        module = []
        for g in genes:
            if x[g].X > 0.5: module.append(g) # well...
        #print sorted(module)
        cov, act_cov, act_sample = 0, 0, 0
        for p in samples:
             cov += int(y[p].X)
             if w[p] > 0:
                 act_sample += 1
                 act_cov += int(y[p].X)



    except GurobiError as e:
        logger.error('Error:', e.message)

    return z, module, e_single, cov, act_cov, act_sample,x,y,m##, p_val_str, p_val_single_str


# Conditional permutation test
def conditional_permutation_test(module, eventToCases, profile, samples, num_permutations):
        # Compute the score for current module
        cases_by_event = [eventToCases[event] for event in module]
        real_W = SuperW(cases_by_event, profile, samples)

        # Find the event with the least significant (highest) P-value
        worst_count, worst_event = 0, None
        for i, (event, cases) in enumerate(zip(module, cases_by_event)):
                # Compute an empirical P-value by shuffling the mutations in the
                # i^th event
                count = 0
                for _ in range(num_permutations):
                        cases_by_event[i] = random.sample(samples, len(cases))
                        count += SuperW(cases_by_event, profile, samples) >= real_W

                # Reset the i^th event to its true cases
                cases_by_event[i] = cases

                # Take only the worst event
                if count > worst_count:
                        worst_count = count
                        worst_event = event

        return worst_event, worst_count

# SuperDendrix weight function
def SuperW(cases_by_event, profile, samples):
        mut_samples = set(s for cases in cases_by_event for s in cases)
        pos = sum(profile[s] for s in mut_samples if profile[s] > 0)
        neg = sum(abs(profile[s]) for cases in cases_by_event for s in cases)
        return 2 * pos - neg

def build_obj(x, y, w, samples_to_genes, args):
    pos_samples = [p for p in samples if w[p] > 0]
    obj_cover_target = quicksum(y[p] * w[p] for p in pos_samples)
    obj_mutex_target = quicksum(quicksum(x[g] * w[p] for g in samples_to_genes[p]) for p in pos_samples) - quicksum(y[p] * w[p] for p in pos_samples)
    obj_neg_muts = quicksum(quicksum(x[g] * -w[p] for g in samples_to_genes[p]) for p in set(samples) - set(pos_samples))
    obj = obj_cover_target
    obj -= obj_mutex_target
    obj -= obj_neg_muts
    return obj

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
