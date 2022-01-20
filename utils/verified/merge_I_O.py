# merge inactivating and other mutations for a gene that co-occur in a cell line
# we keep inactivating and remove other if they co-occur

import copy, sys
from i_o import load_mutation_data

input_features = sys.argv[1]
output_file = sys.argv[2]

print(input_features)
print(output_file)

m, n, genes, patients, mutToSamples, sampleToMut = load_mutation_data(f)

print(m)
print(n)
new_sampleToMut = copy.deepcopy(sampleToMut)
for s in new_sampleToMut.keys():
    for m in sampleToMut[s]:
        if m.endswith("O_MUT"):
            if (m[:-5]+"I_MUT" in new_sampleToMut[s]):
                new_sampleToMut[s].remove(m)


of = open(output_file, 'w')
of.write('#Sample\tEvents\n')
for i in new_sampleToMut:
    of.write(i + '\t' + '\t'.join(list(new_sampleToMut[i])) + '\n')


