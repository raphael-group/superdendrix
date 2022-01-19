import csv, sys, argparse

# rows: samples, cols: mutations

def get_parser():
    d = "convert sample to event format in SuperDendrix events to binary matrix for REVEALER"
    parser = argparse.ArgumentParser(description=d,fromfile_prefix_chars="@")
    parser.add_argument('-i', '--events', required=True, help='input event file for conversion')
    parser.add_argument('-o', '--output', required=True, help='output_prefix')
    return parser

args = get_parser().parse_args(sys.argv[1:])
sample_to_events = dict()
rows = []
with open(args.events,'r') as tsvin:#, open(str(args.output)+'.gct', 'w') as csvout:
    tsvin = csv.reader(tsvin, delimiter='\t')
    for row in tsvin:
        if row[0]=="#Sample":
            header=row
        else:
            sample_to_events[row[0]] = row[1:]
            rows.append(row)

events = set()
for k,v in sample_to_events.items():
    for val in v:
        events.add(val)
samples = set(sample_to_events.keys())
samples = list(samples)
events = list(events)

print(len(events))
print(len(samples))
outputrows = []

for i in range(len(samples)):
    outputrow = []
    outputrow.append(str(samples[i]))
    for j in range(len(events)):
        x = 1 if events[j] in sample_to_events[samples[i]] else 0
        outputrow.append(str(x))
    outputrows.append(outputrow)

#for i in range(len(events)):
#    outputrow = []
#    outputrow.append(str(events[i]))
#    #outputrow.append(str(events[i]))
#    for j in range(len(samples)):
#        x = 1 if events[i] in sample_to_events[samples[j]] else 0
#        outputrow.append(str(x))
#    outputrows.append(outputrow)

for i in range(len(samples)):
    if samples[i][0].isdigit():
        samples[i] = "X" + samples[i]

of = open(args.output+".tsv",'w')
#of.write("#1.2\n")
#of.write(str(len(events)) + "\t" + str(len(samples)) + "\n")
#of.write("Name\tDescription\t" + "\t".join(samples)+"\n")
of.write("Sample\t" + "\t".join(events))
for k in range(len(outputrows)):
    of.write("\n" + "\t".join(outputrows[k]))
of.close()
