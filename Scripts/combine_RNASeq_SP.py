import argparse
import glob
import numpy as np

###################################################################################
parser = argparse.ArgumentParser(description= 'Combine RNA-Seq bedgraph files from S.pneumoniea' + '\n'
								'Usage:' + '\t' + 'position.py <options> -i -o')

#required files:
parser.add_argument('-i', dest='inpath', help='path to RNA-Seq bedgraph files', required=True)
parser.add_argument('-o', dest='outpath', help='outpath', required=True)

args = parser.parse_args()

inpath = args.inpath
outpath = args.outpath
outfile = outpath + 'combined_316_RNASeqND0min_Streptococcus_pneumoniae_genomeCoverage.bedgraph'

##################################################################################



bedgraphFiles = glob.glob(inpath + '*316*.bedgraph')
# print bedgraphFiles

bedgraphFile1 = ''
bedgraphFile2 = ''
bedgraphFile3 = ''
bedgraphFile4 = ''

for i in range(len(bedgraphFiles)):
	locals()["bedgraphFile"+str(i+1)] = ''.join(bedgraphFiles[i])


uniqueReads = 0
chrom = 'NC_003028.3'

with open(bedgraphFile1, 'r') as f1, open(bedgraphFile2, 'r') as f2,\
	open(bedgraphFile3, 'r') as f3, open(bedgraphFile4, 'r') as f4:

	ur = []
	uniqueReads1 = int(f1.readline().split(' ')[2].split('"')[1])
	uniqueReads2 = int(f2.readline().split(' ')[2].split('"')[1])
	uniqueReads3 = int(f3.readline().split(' ')[2].split('"')[1])
	uniqueReads4 = int(f4.readline().split(' ')[2].split('"')[1])

	ur.append([uniqueReads1,uniqueReads2,uniqueReads3,uniqueReads4])
	uniqueReads = int(np.rint(np.mean(ur[0])))
	print ur
	print uniqueReads

	with open(outfile, 'w') as o:
		o.write("track type=bedGraph name=\"" + str(uniqueReads) + "\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n")


		line1 = f1.readline()
		line2 = f2.readline()
		line3 = f3.readline()
		line4 = f4.readline()

		while zip(line1,line2,line3,line4):

			position = int(line1.split()[1])

			coverageCombined = int(np.rint(np.mean([int(line1.split()[2]),int(line2.split()[2]),int(line3.split()[2]),int(line4.split()[2])])))



			o.write(chrom + "\t" + str(position) +  "\t" + str(coverageCombined) + "\n" )

			line1 = f1.readline()
			line2 = f2.readline()
			line3 = f3.readline()
			line4 = f4.readline()
