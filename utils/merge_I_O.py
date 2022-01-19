# merge inactivating and other mutations for a gene that co-occur in a cell line
# we keep inactivating and remove other if they co-occur


import copy, sys
from i_o import load_mutation_data

f = sys.argv[1]
output = sys.argv[2]

#f = "../data/20Q2/features/CERES_oncoKB_MUT.tsv"
#f = "../data/20Q2/features/CERES_oncoKB_CT.tsv"
#f = "../data/20Q2/features/CERES_full_MUT.tsv"
#f = "../data/20Q2/features/CERES_full_CT.tsv"


#f = "../data/Sanger_19Q1/features/SANGER_full_CT.tsv"
##f = "../data/Sanger_19Q1/features/SANGER_full_MUT.tsv"
##f = "../data/Sanger_19Q1/features/SANGER_oncoKB_CT.tsv"
##f = "../data/Sanger_19Q1/features/SANGER_oncoKB_MUT.tsv"

#output = "../data/20Q2/features/IOmerged/CERES_oncoKB_MUT_merged.tsv"
#output = "../data/20Q2/features/IOmerged/CERES_oncoKB_CT_merged.tsv"
#output = "../data/20Q2/features/IOmerged/CERES_full_MUT_merged.tsv"
#output = "../data/20Q2/features/IOmerged/CERES_full_CT_merged.tsv"

#output = "../data/Sanger_19Q1/features/IOmerged/SANGER_full_CT_merged.tsv"
##output = "../data/Sanger_19Q1/features/IOmerged/SANGER_full_MUT_merged.tsv"
##output = "../data/Sanger_19Q1/features/IOmerged/SANGER_oncoKB_CT_merged.tsv"
##output = "../data/Sanger_19Q1/features/IOmerged/SANGER_oncoKB_MUT_merged.tsv"

print(f)
print(output)

m, n, genes, patients, mutToSamples, sampleToMut = load_mutation_data(f)

print(m)
print(n)
new_sampleToMut = copy.deepcopy(sampleToMut)
for s in new_sampleToMut.keys():
    for m in sampleToMut[s]:
        if m.endswith("O_MUT"):
            if (m[:-5]+"I_MUT" in new_sampleToMut[s]):
                new_sampleToMut[s].remove(m)


of = open(output, 'w')
of.write('#Sample\tEvents\n')
for i in new_sampleToMut:
    of.write(i + '\t' + '\t'.join(list(new_sampleToMut[i])) + '\n')


