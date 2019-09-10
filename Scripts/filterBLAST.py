# remove all predicted Terminators with bitscores higher than 30

import argparse

#######################################################################
#######################################################################
parser = argparse.ArgumentParser(description= 'Look for predictions up to 150 nucleotides downstream of genes' + '\n'
								'Usage:' + '\t' + 'position.py <options> -pos -neg -all -gene -term -o')

#required files:
parser.add_argument('-term', dest='predictedTerminators', help='input predicted Terminators', required=True)
parser.add_argument('-blast', dest='blastFile', help='input BLAST', required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)

args = parser.parse_args()

predictedTerminators = args.predictedTerminators
blastFile = args.blastFile
outpath = args.outpath

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
		out.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(item.split('_')[0]) + '\t' + str(int(item.split('_')[0])+1) + '\t' + str(item) +'\n')

