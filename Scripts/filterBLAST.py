# removes all predicted Terminators with bitscores higher than 30
# writes file of predicted terminators of certain length (default 100)
# writes files for embedding the terminator sequences in (500 before and after predicted terminator)

import argparse
import os.path
import math

#######################################################################
#######################################################################
# methods for checking parsed file types

def checkBedFormat(v):
	b = os.path.splitext(v)[1][1:].lower()
	if b != 'bed':
		raise argparse.ArgumentTypeError('bed format file type expected')
	else:
		return v
		
def checkInt(v):
	v = int(v)
	if v < 0:
		raise argparse.ArgumentTypeError('positive Integer value expected')
	if isinstance(v, int):
		return v
	else:
		raise argparse.ArgumentTypeError('Integer value expected')

#######################################################################
#######################################################################

parser = argparse.ArgumentParser(description= 'Filter predicted positives with BLAST bitscores over 30' + '\n'
								'Usage:' + '\t' + 'filterBLAST.py <options> -term -blast -o' +'\n'
								'optional:' + '\t' + '-l')

#required files:
parser.add_argument('-term', dest='predictedTerminators', help='input predicted Terminators directory', type=checkBedFormat, required=True)
parser.add_argument('-blast', dest='blastFile', help='input BLAST directory', required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)
#optional
parser.add_argument('-l', dest='lengthTerminator', help='length of artificial terminator, default:100', type=checkInt, nargs='?', default=100)

args = parser.parse_args()

predictedTerminators = args.predictedTerminators
blastFile = args.blastFile
outpath = args.outpath
lengthTerminator = args.lengthTerminator

l1 = lengthTerminator * 0.16666666666666664
l2 = lengthTerminator - l1

l1 =  int(math.ceil(l1))
l2 = int(math.floor(l2))

organism = ''
chrom = ''
chrom2 = ''
plasmid = ''
lengthGenome = 0


if "BS" in predictedTerminators:
	organism = 'B.subtilis'
	chrom = 'NC_000964.3/1-4215606'
	lengthGenome = 4215606
	brev = 'BS'
if "EF" in predictedTerminators:
	organism = 'E.faecalis'
	if 'chrom' in predictedTerminators:
		chrom = "NC_004668.1"
		plasmid = 'Chromosome'
		lengthGenome = 3218031
		brev = 'EF_chrom'
	if 'pl1' in predictedTerminators:
		chrom = "NC_004669.1"
		plasmid = 'Plasmid1'
		lengthGenome = 66320
		brev = 'EF_pl1'
	if 'pl2' in predictedTerminators:
		chrom = "NC_004671.1"
		plasmid = 'Plasmid2'
		lengthGenome = 57660
		brev = 'EF_pl2'
	if 'pl3' in predictedTerminators:
		chrom = "NC_004670.1"
		plasmid = 'Plasmid3'
		lengthGenome = 17963
		brev = 'EF_pl3'
if "LM" in predictedTerminators:
	organism = 'L.monocytogenes'
	chrom = "NC_003210.1"
	lengthGenome = 2944528
	brev = 'LM'
if 'SP' in predictedTerminators:
	organism = 'S.pneumoniae'
	chrom = "NC_003028.3"
	chrom2 = chrom
	lengthGenome = 2160842
	brev = 'SP'
	
print('\n' + str(organism) + ' ' + str(plasmid))


outfiles = [outpath + brev + '_predTerm_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.bed',
			outpath + brev + '_500front_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.bed',
			outpath + brev + '_500back_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.bed']


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
print("Predicted Terminators before: " + str(len(distanceCoords)))

print("With Bitscores over 30: " + str(len(bitscoreOver30Coords)))

# difference of two sets (A-B): elements only in A but not in B
withoutOver30 = (distanceCoords - bitscoreOver30Coords)

print("Predicted Terminators Without Bitscores over 30: " + str(len(withoutOver30)))


outfiles1 = [open(i, 'w') for i in outfiles]

for item in withoutOver30:
	coord = int(item.split()[-1].split('_')[0])
	strand = item.split()[-1].split('_')[-1]

	if coord-l2+500 > 0 and coord+l2+500 <= lengthGenome:
		if strand == '-':
			outfiles1[0].write(str(chrom) + '\t' + str(coord-l2) + '\t' + str(coord+l1) + '\t' + str(item) +'\n')

			outfiles1[1].write(str(chrom) + '\t' + str(coord-l2-1-500) + '\t' + str(coord-1-l2) + '\t' + str(item) +'\n')
			outfiles1[2].write(str(chrom) + '\t' + str(coord+l1+1) + '\t' + str(coord+l1+1+500) + '\t' + str(item) +'\n')
		else:
			outfiles1[0].write(str(chrom) + '\t' + str(coord-l1) + '\t' + str(coord+l2) + '\t' + str(item) +'\n')

			outfiles1[1].write(str(chrom) + '\t' + str(coord-l1-1-500) + '\t' + str(coord-1-l1) + '\t' + str(item) +'\n')
			outfiles1[2].write(str(chrom) + '\t' + str(coord+l2+1) + '\t' + str(coord+l2+1+500) + '\t' + str(item) +'\n')