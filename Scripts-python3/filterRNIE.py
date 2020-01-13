# filters RNIE output by score, write file with RNIE scores over 20 and make scatterplots

import argparse
import os.path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import collections
import sys
import os.path

#######################################################################
#######################################################################

def checkGffFormat(v):
	b = os.path.splitext(v)[1][1:].lower()
	if b != 'gff':
		raise argparse.ArgumentTypeError('gff format file type expected')
	else:
		return v

#######################################################################
#######################################################################
parser = argparse.ArgumentParser(description= 'Filter RNIE output by score (>20.0)' + '\n'
								'Usage:' + '\t' + 'filterRNIE.py <options> -i -o')

#required files:
parser.add_argument('-i', dest='rnieFile', help='path to RNIE output in gff format', type=checkGffFormat, required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)


args = parser.parse_args()

rnieFile = args.rnieFile
outpath = args.outpath

outfile = outpath + 'filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed'
plot = outpath + 'filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.png'

#######################################################################
#######################################################################

RNIEdict = {}
xs = []
ys = []
RNIEover20X = []
RNIEover20Y = []

with open(rnieFile, 'r') as rnie:

	for line in rnie:
		coord = int(line.split('_')[0])
		x = float(line.split('_')[1])
		y = float(line.split('_')[2])
		strand = line.split('_')[3].split()[0]
		score = float(line.split()[5])
		xs.append(x)
		ys.append(y)


		if coord not in RNIEdict and score > 20.0:
			RNIEdict[coord] = [x,y,strand,score]
			RNIEover20X.append(x)
			RNIEover20Y.append(y)


with open(outfile, 'w') as out:
	for value, key in list(RNIEdict.items()):
		out.write(str(value) + '\t' + str(key[0]) + '\t' + str(key[1]) + '\t' + str(key[2]) + '\t' + str(key[3]) + '\n')


npX = np.array(xs)
npY = np.array(ys)
npRNIEover20X = np.array(RNIEover20X)
npRNIEover20Y = np.array(RNIEover20Y)

print('RNIE scores over 20.0: ' + str(len(RNIEover20X)))

#######################################################################
#######################################################################
fig, ax = plt.subplots(1,1, figsize=(6,6), dpi=120)
plt.title('Filter RNIE scores', fontsize=14)

plt.scatter(npY+1, npX+1, s=15)
plt.scatter(npRNIEover20Y+1, npRNIEover20X+1, label='RNIE over 20.0', c='red', marker='x', s=12 )
plt.xscale('log')
plt.yscale('log')
plt.xlim(1, 100000)
plt.ylim(1, 100000)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Avg. RNA-Seq count', fontsize=14)
plt.ylabel('Max. Term-Seq count', fontsize=14)
plt.grid(True)
ax.set_axisbelow(True)

plt.legend(prop={'size': 14})

plt.savefig(plot,dpi=300)
