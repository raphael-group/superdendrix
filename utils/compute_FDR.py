import sys, time, argparse, statistics, numpy
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

stats= importr('stats')

def get_parser():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@')
    parser.add_argument('-i', '--input_file', required=True, help='input concat SD results')
    parser.add_argument('-p', '--pval_index', required=False, default=4, type=int, help='pvalue index in input file')
    parser.add_argument('-o', '--output_file', required=True, help='name of output file')
    return parser

args = get_parser().parse_args(sys.argv[1:])
header=[]
rows = []
seenheader = False
pi = args.pval_index
ps = []
with open(args.input_file) as f:
    for r in f:
        if not seenheader:
            row = r.strip("\n").split("\t")
            row.append('FDR')
            header=row
            seenheader = True
        else:
            row = r.strip("\n").split("\t")
            rows.append(row)
            if row[pi] == "na":
                ps.append(float(1.0))
            else:
                ps.append(float(row[pi]))

p_adjust = list(stats.p_adjust(FloatVector(ps), method='BH'))

of = open(args.output_file, 'w')
print("header: ",header)
print("check that this is pval: ",header[pi])
print("first result row: ",rows[0])
of.write('\t'.join(header))
for i in range(len(rows)):
    of.write('\n' + '\t'.join(rows[i]) + '\t' + str(p_adjust[i]))
of.close()
