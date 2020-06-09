#!/usr/bin/env python

# Load required modules
# Update: GMM 0.20



# 0516:  do not exclude top/bottom : original

# 0530: ratio based on comp size (background vs significant)

# 1125: ratio based on means (higher) / (lower)
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
parser.add_argument('-ef', '--events_file', type=str, required=True)
parser.add_argument('-o', '--output_file', type=str, required=True)
parser.add_argument('-nc', '--num_cores', type=int, default=1, required=False)
parser.add_argument('--output_models', action='store_true', default=False, required=False,
                    help='Write one file per GMM')
parser.add_argument('-cv', '--covariance_type', type=str, required=False, default='spherical',
                    choices=['tied', 'diag', 'full', 'spherical'])
parser.add_argument('-rs', '--rand_seed', type=int, default=1764, required=False)
args = parser.parse_args(sys.argv[1:])

np.random.seed(args.rand_seed)
# Load the events file
eventToCases, mutation_samples = load_events(args.events_file, verbose=1)

# Load the profile file
profiles, profile_samples = load_profiles(args.profile_file, verbose=1)
sampleToIndex = dict(zip(profile_samples, range(len(profile_samples))))
events = sorted(profiles.keys())
profile_sample_set = set(profile_samples)
sampleToIndex = dict(zip(profile_samples, range(len(profile_samples))))


def event_freq(name, ty):
    event = name.split('_')[0] + '_' + ty.upper()
    if event in eventToCases:
        LL = sum(weighted_profiles[name][sampleToIndex[s]]
                 for s in eventToCases[event] if s in profile_sample_set)
        num = len(eventToCases[event] & profile_sample_set)
    else:
        LL = '--'
        num = '--'
    return num, LL

# Fit GMM function

def fit_gmm(profile):
    # Featurize the data
#    mean = np.mean(profile)
#    std = np.std(profile)
#    zscores = (profile - mean)/std
##    outlier = [abs(z) > 7 for z in zscores]  # zscore
##    Max = max(profile)
##    Min = min(profile)
##    outlier = [(p >= Max) or (p <= Min) for p in profile]  # exclude the largest and the smallest
##    condition = [(x or y) for x,y in zip(np.isnan(profile),outlier)]

    nan_indices = np.where(np.isnan(profile))[0]
##    excluded_indices = np.where(condition)[0]

##    X_both = np.array([[profile[p]] for p in range(len(profile)) if not condition[p]])
##    X_nan = np.array([[p] for p in profile if not np.isnan(p)])
    X = np.array([[p] for p in profile if not np.isnan(p)])


    # Choose number of components
    lowest_bic = np.infty
    best_ll = np.infty
    bic = []
    ll = []
    n_components_range = range(1, 3)
    for K in n_components_range:
        # Fit a Gaussian mixture with EM
        #gmm = GMM(n_components=K,covariance_type='spherical',)
        gmm = GMM(n_components=K, random_state= args.rand_seed, n_init=20, max_iter=500, covariance_type=args.covariance_type)
        gmm.fit(X)
        bic.append(gmm.bic(X))
        ll.append(gmm.score(X))
#        if bic[-1] < lowest_bic:
#            lowest_bic = bic[-1]
#            best_ll = ll[-1]
#            best_gmm = gmm
#            if len(bic) > 1:
#                higher_bic = bic[0]
#                best_gmm.BIC2 = higher_bic
    best_gmm = gmm
    best_gmm.BIC1 = bic[0]
    best_gmm.BIC2 = bic[1]
    # Store X so we can use it later
    best_gmm.X = X
    best_gmm.nan_indices = nan_indices
##    best_gmm.excluded_indices = excluded_indices
    best_gmm.LL1 = ll[0]
    best_gmm.LL2 = ll[1]
##    best_gmm.excluded_num = sum(condition)

    return best_gmm

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
models = dict()
excluded = set()
rows = []
for name, gmm in zip(events, gmms):
    # Compute new weights as posterior probability of
    # each sample belonging to the GMM component with highest
    # mean
    means = gmm.means_[:, 0]
#        var = gmm.covariances_[:, 0]
    var = gmm.covariances_ # for spherical 
    resp = gmm.predict_proba(gmm.X)
    labels = gmm.predict(gmm.X)
    nnans = len(gmm.nan_indices)
    weights = gmm.weights_
#        if nnans > 0:
#            print(name, nnans)
'''
    if np.max(weights) == weights[0]:
        k, l = 0, 1 # k is higher
    else:
        k, l = 1, 0
'''

    if np.max(means) == means[0]:
        i, j = 0, 1 # i is higher
    else:
        i, j = 1, 0

#        print("###################HERE##################")
#        print("###################HERE##################")
#        print("lower ",resp[:,i])
#        print("higher ",resp[:,j])
#        print("resp ",resp)
#        print("means ",means)
    ## higher / lower
    ##weighted_profiles[name] = list(np.log(resp[:, i] / resp[:, j]))
    
    ##component with higher mean / lower mean
    #rat = np.log(resp[:, l] / resp[:, k])
    rat = np.log(resp[:, i] / resp[:, j])
    toolarge = 0
    toosmall = 0
    for k in range(len(rat)):
        if rat[k] > 20:
            rat[k] = 20
            toolarge += 1
        elif rat[k] < -20:
            rat[k] = -20
            toosmall += 1

    #weighted_profiles[name] = list(np.log(resp[:, j] / resp[:, i]))
    weighted_profiles[name] = list(rat)
    # Add back the NaNs as having weight zero
    for idx in reversed(gmm.nan_indices):
        weighted_profiles[name].insert(idx, 0)

    # Sum the LLRs for samples with mutations/amplifications/deletions
    # in the profile gene
    ##numMut, mutLL = event_freq(name, 'MUT')
    ##numAmp, ampLL = event_freq(name, 'AMP')
    ##numDel, delLL = event_freq(name, 'DEL')

    # Record the row and the model
#        m1, v1, n1 = means[j], var[j][0], sum(labels == j) # comp with lower mean
    
#        m2, v2, n2 = means[i], var[i][0], sum(labels == i) # comp with higher mean

    m1, v1, n1 = means[j], var[j], sum(labels == j) # comp with lower mean
    
    m2, v2, n2 = means[i], var[i], sum(labels == i) # comp with higher mean

    (w1, w2) = (gmm.weights_[j],gmm.weights_[i])
    
    rows.append((name, m1, v1, n1,w1, m2, v2, n2,w2, nnans, str(gmm.BIC1),str(gmm.BIC2), str(gmm.BIC1 - gmm.BIC2), str(toolarge), str(toosmall), str(gmm.LL2)))##,
##        rows.append((name, m1, v1, n1, m2, v2, n2, nnans, str(gmm.BIC), str(gmm.BIC2 - gmm.BIC),str(gmm.excluded_num), str(toolarge), str(toosmall)))##,
                 ##numMut, mutLL, numAmp, ampLL, numDel, delLL))

    models[name] = gmm
else:
    excluded.add(name)

print('\t- %s profiles with >= 1 components' % len(weighted_profiles))
print('\t- %s profiles excluded' % len(excluded))

# Output the models and then the new profiles to file
output_prefix = os.path.splitext(args.output_file)[0]




with open(output_prefix + '-profile-summary.tsv', 'w') as OUT:
    OUT.write(
        '#Profile name\tMean 1\tVar 1\tN1\tWeight 1\tMean 2\tVar 2\tN2\tWeight 2\tNans\tk1BIC\tk2BIC\tdiff(BIC)\tnum_+capped\tnum_-capped\tk2LL\n')
    OUT.write('\n'.join('\t'.join(map(str, r)) for r in rows))


with open(output_prefix + '-events-w-two-components.txt', 'w') as OUT:
    OUT.write('#Events with two-components in the mixture model\n')
    OUT.write('\n'.join(sorted(models.keys())))

if args.output_models:
    params = dict(vars(args))
    output = dict(models=models, params=params)
    joblib.dump(output, output_prefix + '-models.pkl')

profile_names = sorted(weighted_profiles.keys())
profile_matrix = [[profile_samples[i]] + [weighted_profiles[n][i]
                                          for n in profile_names] for i in range(len(profile_samples))]

with open(args.output_file, 'w') as OUT:
    OUT.write('#\t%s\n' % '\t'.join(profile_names))
    OUT.write('\n'.join('\t'.join(map(str, r)) for r in profile_matrix))
