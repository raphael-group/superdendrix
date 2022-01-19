# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 19:35:37 2020

@author: pty01
"""
import numpy as np
import copy

d_profile_summary = "../data/21Q1/profiles/GMM/CERES_zscores_2C_with_resp-profile-summary.tsv"
d_twoC = "../data/21Q1/profiles/GMM/CERES_zscores_2C_with_resp_2Cscores.tsv"
d_rawscore = "../data/21Q1/profiles/CERES_zscores.tsv"


profile_summary = dict()
with open(d_profile_summary, 'r') as f:
    arrs = [l.rstrip('\n').split('\t')
            for l in f]
    header = arrs.pop(0)[1:]
    print(header)
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

#np.nanmax([np.nan,1,2,3])
#np.nanmax(list(twoC['KCNK13 (56659)'].values()))

genes = list(profile_summary.keys())
samples = list(twoC[genes[0]].keys())
print(genes[0])
print(samples)

new_twoC = copy.deepcopy(twoC)


for g in range(len(profile_summary)):
    gene = genes[g]
    c = 0
    
    maxscore = np.nanmax(list(twoC[gene].values()))
    minscore = np.nanmin(list(twoC[gene].values()))
    
    for s in range(len(samples)):
        sample = samples[s]
        
        if not np.isnan(rawscore[gene][sample]):
            if (rawscore[gene][sample] > float(profile_summary[gene]["Mean 2"])) and (twoC[gene][sample] < 0):
                new_twoC[gene][sample] = maxscore
                c += 1
            if (rawscore[gene][sample] < float(profile_summary[gene]["Mean 1"])) and (twoC[gene][sample] > 0):
                new_twoC[gene][sample] = minscore
                c += 1
    if g % 1000 == 0:
        print(gene)
        print(c)
        print(g)
    profile_summary[gene]["crossNum"] = c

#profile_summary['BRAF (673)']


outputfile = "../data/20Q2/profiles/GMM/corrected_CERES_zscores_2C_with_resp-profile-summary.tsv"
of = open(outputfile, 'w')
of.write("profile" + "\t" + "\t".join(list(list(profile_summary.values())[0].keys())) + "\n")
for k,v in profile_summary.items():
    of.write(k + "\t")
    of.write("\t".join(list(map(str,list(v.values())))) + "\n")
of.close()    



outputfile2 = "../data/20Q2/profiles/GMM/corrected_CERES_zscores_2C_with_resp_2Cscores.tsv"
of2 = open(outputfile2, 'w')

of2.write("Sample" + "\t" + "\t".join(list(new_twoC.keys())) + "\n")

for sample in samples:
    curr = [sample]
    for gene in list(new_twoC.keys()):
        curr.append(str(new_twoC[gene][sample]))
    of2.write("\t".join(curr) + "\n")
of2.close()
   



