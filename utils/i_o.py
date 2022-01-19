#!/usr/bin/env python

# Load required modules
import sys, os, logging
from collections import defaultdict
import numpy as np
# Set up a logger
FORMAT = '%(asctime)s %(filename)-15s %(levelname)-10s: %(message)s'
logging.basicConfig(format=FORMAT)

def getLogger(verbosity=logging.INFO):
    logger = logging.getLogger(__name__)
    logger.setLevel(verbosity)
    return logger

###############################################################################
# Functions for loading mutation data
def load_mutation_data(filename, patients=None, geneFile=None, mut_only = False, min_freq = 0, max_freq = float("inf"), verbose = False):
    """Loads the mutation data in the given file.

    :type filename: string
    :param filename: Location of mutation data file.
    :type patient_file: string
    :param patient_file: Location of patient (whitelist) file.
    :type gene_file: string
    :param gene_file: Location of gene (whitelist) file.
    :rtype: Tuple

    **Returns:**
      * **m** (*int*) - number of patients.
      * **n** (*int*) - number of genes.
      * **genes** (*list*) - genes in the mutation data.
      * **patients** (*list*) - patients in the mutation data.
      * **geneToCases** (*dictionary*) - mapping of genes to the patients in which they are mutated.
      * **patientToGenes** (*dictionary*) - mapping of patients to the genes they have mutated.
    """
    # # Load the whitelists (if applicable)
    # if patientFile:
    #     with open(patientFile) as f:
    #         patients = set( l.rstrip().split()[0] for l in f if not l.startswith("#") )
    # else:
    #     patients = None

    if geneFile:
        with open(geneFile) as f:
            genes = set( l.rstrip().split()[0] for l in f if not l.startswith("#") )
    else:
        genes = set()

    # Parse the mutation matrix
    from collections import defaultdict
    geneToCases, patientToGenes = defaultdict(set), defaultdict(set)
    with open(filename) as f:
        arrs = [ l.rstrip().split("\t") for l in f if not l.startswith("#") ]
        for arr in arrs:
            patient, mutations = arr[0], set(arr[1:])
            if not patients or patient in patients:
                if genes: mutations &= genes
                #else: genes |= mutations

                patientToGenes[patient] = mutations
                for gene in mutations:
                    geneToCases[gene].add(patient)


    # Format and return output
    genes, patients = list(geneToCases.keys()), list(patientToGenes.keys())

    if mut_only:
        toRemove = [ g for g in genes if not g.endswith("MUT") ]
        if verbose: print("Warning:", len(toRemove), "events removed because no mutations")
        for g in toRemove:
            for p in geneToCases[g]:
                patientToGenes[p].remove(g)
            del geneToCases[g]
            genes.remove(g)

    toRemove = [ g for g in genes if len(geneToCases[g]) < min_freq or len(geneToCases[g]) > max_freq ]
    if verbose: print("Warning:", len(toRemove), "events removed because of too low/high frequency")
    for g in toRemove:
        for p in geneToCases[g]:
            patientToGenes[p].remove(g)
        del geneToCases[g]
        genes.remove(g)

    m, n = len(genes), len(patients)

    return m, n, genes, patients, geneToCases, patientToGenes

 # Load event file
def load_events(mutation_file, verbose=0):
    if verbose > 0:
        print('* Loading mutations...')
    with open(mutation_file, 'r') as IN:
        arrs = [l.rstrip('\n').split('\t')
                for l in IN if not l.startswith('#')]
        eventToCases, mutation_samples = defaultdict(set), set()
        for arr in arrs:
            sample = arr[0]
            mutation_samples.add(sample)
            for event in arr[1:]:
                eventToCases[event].add(sample)
    if verbose > 0:
        print('\t- Events: %s' % len(eventToCases))
    return eventToCases, mutation_samples


# Load the profile file
def load_profiles(profile_file, sample_whitelist=None, verbose=0):
    if verbose > 0:
        print('* Loading profiles...')
    with open(profile_file, 'r') as IN:
        if profile_file.endswith(".tsv"):
            delim = "\t"
        elif profile_file.endswith(".csv"):
            delim = ","
        arrs = [l.rstrip('\n').split(delim)
                for l in IN if not l.startswith('#')]
#        arrs = [l.rstrip('\n').split(',')
#                for l in IN]
        print("#####################")
        profile_events = arrs.pop(0)[1:]
        profile_samples = [arr[0] for arr in arrs]
        #print(profile_samples)
        #print(profile_events)
        profile_matrix = np.array(
            [[float(a) if (a.lower() != 'na') and (a != '') else np.nan for a in arr[1:]] for arr in arrs])

        # Restrict to given list of samples in the whitelist
        if sample_whitelist is not None:
            sample_indices = [i for i, s in enumerate(
                profile_samples) if s in sample_whitelist]
            samples = [profile_samples[i] for i in sample_indices]
            profiles = dict((e, profile_matrix[sample_indices, i])
                            for i, e in enumerate(profile_events))
        else:
            samples = profile_samples
            profiles = dict((e, profile_matrix[:, i])
                            for i, e in enumerate(profile_events))

    if verbose > 0:
        print('\t- Number of profiles:', len(profiles))
        print('\t- Samples in common:', len(samples))

    return profiles, samples

# Load one differential dependency
def load_single_profile(profile_fn, gene):
    # Generate/load the target profile
    w = {} # weights
    global numNanInf
    numNanInf = 0
    with open(profile_fn) as f:
        line = f.readline()
        delim = "\t" if profile_fn.endswith(".tsv") else ","
        line2 = line.rstrip().split(delim)
        line2 = [g.split("_")[0] for g in line2]
        #ind = line.rstrip().split(delim).index(gene)# + 1
        ind = line2.index(gene)# + 1
        arrs = [ l.rstrip().split(delim) for l in f if not l.startswith("#") ]
        #arrs = [ re.findall(r"[-\w']+", l) for l in f if not l.startswith("#") ]
        #for arr in [arr for arr in arrs if arr[0] in patients]:
        for arr in arrs:
            if arr[ind] == "":
                w[arr[0]] = 0.0
            else:
                w[arr[0]] = float(arr[ind])
    return w



 # Load cancer type information for each cell line (CL\tCT)
def load_cancertype(ct_file, verbose=0, IC=True, th=10):
    if verbose > 0:
        print('* Loading cancer types...')
    with open(ct_file, 'r') as IN:
        arrs = [l.rstrip('\n').split('\t')
                for l in IN if not l.startswith('#')]
        CTtoCL = defaultdict(set)
        CLtoCT = defaultdict(set)
        for arr in arrs:
            CL = arr[0]
            if IC and CL[0].isdigit():
                CL = 'X' + CL
            CT = arr[1]

            CTtoCL[CT].add(CL)
            CLtoCT[CL].add(CT)
        ex = 0

    if verbose > 0:
        print('\t- Cancer Types: %s' % len(CTtoCL))
        print('\t- Cancer Types excluded: %s' % ex)
    return CTtoCL,CLtoCT

