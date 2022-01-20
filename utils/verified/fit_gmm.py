#!/usr/bin/env python

# Load required modules
# Update: GMM 0.20

import sys
import os
import argparse
import numpy as np
import itertools
import multiprocessing as mp
from sklearn.mixture import GaussianMixture as GMM
from sklearn.externals import joblib
from collections import defaultdict

import warnings
warnings.filterwarnings("ignore")

# Load our modules
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
from i_o import *

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pf', '--profile_file', type=str, required=True)
parser.add_argument('-o', '--output_directory', type=str, required=True)
parser.add_argument('-nc', '--num_cores', type=int, default=1, required=False)
parser.add_argument('--output_models', action='store_true', default=False, required=False,
                    help='Write one file per GMM')
parser.add_argument('-cv', '--covariance_type', type=str, required=False, default='spherical',
                    choices=['tied', 'diag', 'full', 'spherical'])
parser.add_argument('-rs', '--rand_seed', type=int, default=1764, required=False)

args = parser.parse_args(sys.argv[1:])

# Set random seed
np.random.seed(args.rand_seed)

# Load the profile file
profiles, profile_samples = load_profiles(args.profile_file, verbose=1)
events = sorted(profiles.keys())
sampleToIndex = dict(zip(profile_samples, range(len(profile_samples))))

# Fit GMM function

def fit_gmm(profile):
    # Featurize the data
    nan_indices = np.where(np.isnan(profile))[0]
    X = np.array([[p] for p in profile if not np.isnan(p)])

    # Choose number of components
    lowest_bic = np.infty
    best_ll = np.infty
    bic = []
    ll = []
    n_components_range = range(1, 3)
    for K in n_components_range:
        # Fit a Gaussian mixture with EM
        gmm = GMM(n_components=K, random_state= args.rand_seed, n_init=20, max_iter=500, covariance_type=args.covariance_type)
        gmm.fit(X)
        bic.append(gmm.bic(X))
        ll.append(gmm.score(X))

    gmm_results = gmm
    gmm_results.BIC1 = bic[0]
    gmm_results.BIC2 = bic[1]
    # Store X so we can use it later
    gmm_results.X = X
    gmm_results.nan_indices = nan_indices
    gmm_results.LL1 = ll[0]
    gmm_results.LL2 = ll[1]

    return gmm_results

# Set up the multiprocessing
print('* Fitting GMMs....')
num_cores = args.num_cores if args.num_cores != -1 else mp.cpu_count()
if num_cores != 1:
    pool = mp.Pool(num_cores)
    map_fn = pool.map
else:
    map_fn = map

gmm_args = [profiles[e] for e in events]
gmms = map_fn(fit_gmm, gmm_args)
if num_cores != 1:
    pool.close()
    pool.join()

# Filter the resulting profiles
print('* Filtering GMMs...')
weighted_profiles = dict()
responsibilities = dict()
models = dict()
excluded = set()
rows = []
for name, gmm in zip(events, gmms):
    # Compute new weights as posterior probability of
    # each sample belonging to the GMM component with highest
    # mean
    means = gmm.means_[:, 0]
    var = gmm.covariances_ # for spherical 
    resp = gmm.predict_proba(gmm.X)
    labels = gmm.predict(gmm.X)
    nnans = len(gmm.nan_indices)
    weights = gmm.weights_
#        if nnans > 0:
#            print(name, nnans)

    if np.max(means) == means[0]:
        i, j = 0, 1 # i is higher
    else:
        i, j = 1, 0

    ##component with higher mean / lower mean
    rat = np.log(resp[:, i] / resp[:, j])
    responsibility = resp[:, i]
    toolarge = 0
    toosmall = 0

    # Hard caps on the 2C scores: 20, -20
    for k in range(len(rat)):
        if rat[k] > 20:
            rat[k] = 20
            toolarge += 1
        elif rat[k] < -20:
            rat[k] = -20
            toosmall += 1

    weighted_profiles[name] = list(rat)
    responsibilities[name] = list(responsibility)
    
    # Add back the NaNs as having weight zero
    for idx in reversed(gmm.nan_indices):
        weighted_profiles[name].insert(idx, 0)

    for idx in reversed(gmm.nan_indices):
        responsibilities[name].insert(idx, 0)

    m1, v1, n1 = means[j], var[j], sum(labels == j) # comp with lower mean
    m2, v2, n2 = means[i], var[i], sum(labels == i) # comp with higher mean
    (w1, w2) = (gmm.weights_[j],gmm.weights_[i])
    
    rows.append((name, m1, v1, n1,w1, m2, v2, n2,w2, nnans, str(gmm.BIC1),str(gmm.BIC2), str(gmm.BIC1 - gmm.BIC2), str(toolarge), str(toosmall), str(gmm.LL2)))##,
    models[name] = gmm


print('\t- %s profiles with >= 1 components' % len(weighted_profiles))

# Output the models and then the new profiles to file

with open(os.path.join(args.output_directory, 'gmm_profile_summary.tsv'), 'w') as OUT:
    OUT.write(
        '#Profile_name\tMean_1\tVar_1\tN1\tWeight_1\tMean_2\tVar_2\tN2\tWeight_2\tNans\tk1BIC\tk2BIC\tdiff(BIC)\tnum_+capped\tnum_-capped\tk2LL\n')
    OUT.write('\n'.join('\t'.join(map(str, r)) for r in rows))

profile_names = sorted(weighted_profiles.keys())
profile_matrix = [[profile_samples[i]] + [weighted_profiles[n][i]
                                          for n in profile_names] for i in range(len(profile_samples))]

with open(os.path.join(args.output_directory, '2C_scores.tsv'), 'w') as OUT:
    OUT.write('Sample\t%s\n' % '\t'.join(profile_names))
    OUT.write('\n'.join('\t'.join(map(str, r)) for r in profile_matrix))
