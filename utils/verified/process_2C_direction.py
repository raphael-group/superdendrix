import os, sys
import numpy as np
import copy
import argparse


# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gmm', '--gmm_summary', type=str, required=True)
parser.add_argument('-twoc', '--twoc_scores', type=str, required=True)
parser.add_argument('-raw', '--raw_scores', type=str, required=True)

parser.add_argument('-o', '--output_directory', type=str, required=True)

args = parser.parse_args(sys.argv[1:])

d_profile_summary = args.gmm_summary
d_twoC = args.twoc_scores
d_rawscore = args.raw_scores

summary_file = os.path.join(args.output_directory, 'monotonicity_summary.tsv')
scores_file = os.path.join(args.output_directory, 'monotonic_twoC_scores.tsv')


# Read input files: GMM results, 2C scores, CERES zscores
profile_summary = dict()
with open(d_profile_summary, 'r') as f:
    arrs = [l.rstrip('\n').split('\t')
            for l in f]
    header = arrs.pop(0)[1:]
    for arr in arrs:
        profile_summary[arr[0]] = dict(zip(header, arr[1:]))

twoC = dict()
with open(d_twoC, 'r') as f:
    arrs = [l.rstrip('\n').split('\t')
            for l in f if not l.startswith("#")]
    profile_genes = arrs.pop(0)[1:]
    profile_samples = [arr[0] for arr in arrs]
    
    for i in range(len(profile_genes)):
        twoC[profile_genes[i]] = dict()
        for s in range(len(profile_samples)):
            twoC[profile_genes[i]][profile_samples[s]] = float(arrs[s][1+i])

rawscore = dict()
with open(d_rawscore, 'r') as f:
    arrs = [l.rstrip('\n').split('\t')
            for l in f if not l.startswith("#")]
    profile_genes = arrs.pop(0)[1:]
    profile_samples = [arr[0] for arr in arrs]
    
    for i in range(len(profile_genes)):
        rawscore[profile_genes[i]] = dict()
        for s in range(len(profile_samples)):
            if arrs[s][1+i] == 'NA':
                rawscore[profile_genes[i]][profile_samples[s]] = np.nan
            else:
                rawscore[profile_genes[i]][profile_samples[s]] = float(arrs[s][1+i])


genes = list(profile_summary.keys())
samples = list(twoC[genes[0]].keys())

new_twoC = copy.deepcopy(twoC)
# Ensure that the 2C scores and CERES z-scores are monotonic.
for g in range(len(profile_summary)):
    gene = genes[g]
    c = 0

    maxscore = np.nanmax(list(twoC[gene].values()))
    minscore = np.nanmin(list(twoC[gene].values()))
    
    for s in range(len(samples)):
        sample = samples[s]
        
        if not np.isnan(rawscore[gene][sample]):
            if (rawscore[gene][sample] > float(profile_summary[gene]["Mean_2"])) and (twoC[gene][sample] < 0):
                new_twoC[gene][sample] = maxscore
                c += 1
            if (rawscore[gene][sample] < float(profile_summary[gene]["Mean_1"])) and (twoC[gene][sample] > 0):
                new_twoC[gene][sample] = minscore
                c += 1

    profile_summary[gene]["crossNum"] = c

# Record results
of = open(summary_file, 'w')
of.write("profile" + "\t" + "\t".join(list(list(profile_summary.values())[0].keys())) + "\n")
for k,v in profile_summary.items():
    of.write(k + "\t")
    of.write("\t".join(list(map(str,list(v.values())))) + "\n")
of.close()    

# Save new scores
of2 = open(scores_file, 'w')
of2.write("Sample" + "\t" + "\t".join(list(new_twoC.keys())) + "\n")
for sample in samples:
    curr = [sample]
    for gene in list(new_twoC.keys()):
        curr.append(str(new_twoC[gene][sample]))
    of2.write("\t".join(curr) + "\n")
of2.close()
   



