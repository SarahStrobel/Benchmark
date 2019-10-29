# remove all predicted Terminators with bitscores higher than 30

import argparse
import os.path
#######################################################################
#######################################################################
# methods for checking parsed file types

def checkBedFormat(v):
	b = os.path.splitext(v)[1][1:].lower()
	if b != 'bed':
		raise argparse.ArgumentTypeError('bed format file type expected')
	else:
		return v
#######################################################################
#######################################################################
parser = argparse.ArgumentParser(description= 'Filter predicted positives with BLAST bitscores over 30' + '\n'
								'Usage:' + '\t' + 'filterBLAST.py <options> -term -blast -o')

#required files:
parser.add_argument('-term', dest='predictedTerminators', help='input predicted Terminators', type=checkBedFormat, required=True)
parser.add_argument('-blast', dest='blastFile', help='input BLAST', required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)

args = parser.parse_args()

predictedTerminators = args.predictedTerminators
blastFile = args.blastFile
outpath = args.outpath

organism = ''
chrom = ''
chrom2 = ''
plasmid = ''

if "BS" in predictedTerminators:
	organism = 'B.subtilis'
	chrom = 'NC_000964.3/1-4215606'
	chrom2 = 'gi|255767013|ref|NC_000964.3|'
if "EF" in predictedTerminators:
	organism = 'E.faecalis'
	if 'chrom' in predictedTerminators:
		chrom = "NC_004668.1"
		chrom2 = chrom
		plasmid = 'Chromosome'
	if 'pl1' in predictedTerminators:
		chrom = "NC_004669.1"
		chrom2 = chrom
		plasmid = 'Plasmid1'
	if 'pl2' in predictedTerminators:
		chrom = "NC_004671.1"
		chrom2 = chrom
		plasmid = 'Plasmid2'
	if 'pl3' in predictedTerminators:
		chrom = "NC_004670.1"
		chrom2 = chrom
		plasmid = 'Plasmid3'
if "LM" in predictedTerminators:
	organism = 'L.monocytogenes'
	chrom = "NC_003210.1"
	chrom2 = chrom
# if 'SP' in predictedTerminators:
# 	organism = 'S.pneumoniae'
# 	chrom = "NC_003028.3"
# 	chrom2 = chrom

print '\n' + str(organism) + ' ' + str(plasmid)

outfile = outpath + 'BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.bed'

#######################################################################
#######################################################################

bitscoreOver30Coords = set()
distanceCoords = set()

with open(predictedTerminators,'r') as term, open(blastFile,'r') as blast:

	for line1 in blast:
		if float(line1.split()[11]) > 30.0:
			coordBLAST = str(line1.split()[0])
			bitscoreOver30Coords.add(coordBLAST)

	for line2 in term:
		coordDistance = line2.split()[3]		
		distanceCoords.add(coordDistance)


# print sorted(distanceCoords)
print "Predicted Terminators before: " + str(len(distanceCoords))

print "With Bitscores over 30: " + str(len(bitscoreOver30Coords))

# difference of two sets (A-B): elements only in A but not in B
withoutOver30 = (distanceCoords - bitscoreOver30Coords)

print "Predicted Terminators Without Bitscores over 30: " + str(len(withoutOver30))

with open(outfile, 'w') as out:
	for item in withoutOver30:
		coord = int(item.split()[-1].split('_')[0])
		coord2 = coord + 1
		# print coord
		out.write(str(chrom) + '\t' + str(coord) + '\t' + str(coord2) + '\t' + str(item) +'\n')

