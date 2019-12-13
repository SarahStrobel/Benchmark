# takes output from ratioTermRNASeq.py and makes scatterplots of training data +
# takes output from classification.py (bed files) and adds scatterplots of  
# predicted terminators and predicted negatives + colors according to their distances to annotated gene


import numpy as np
import math
import matplotlib.pyplot as plt
from bisect import bisect_left
from matplotlib.colors import ListedColormap
import argparse
import os.path

#######################################################################
#######################################################################
# methods for checking parsed file types

def checkGffFormat(v):
	b = os.path.splitext(v)[1][1:].lower()
	if not (b == 'gff' or b == 'gff3'):
		raise argparse.ArgumentTypeError('gff or gff3 format file type expected')
	else:
		return v

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
# https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
# requires myList to be sorted (gene annotation should be!) returns closest value to myNumber

def closestGene(myList, myNumber): 
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
       return after
    else:
       return before


#######################################################################
#######################################################################
parser = argparse.ArgumentParser(description= 'Look for predictions up to 150 nucleotides downstream of genes' + '\n'
								'Usage:' + '\t' + 'position.py <options> -pos -neg -all -gff -term -o' + '\n'
								'optional:' + '\t' + '-l')

#required files:
parser.add_argument('-pos', dest='predictedTerminators', help='input predicted Terminators', type=checkBedFormat, required=True)
parser.add_argument('-neg', dest='predictedNegatives', help='input predicted Negatives', type=checkBedFormat, required=True)
parser.add_argument('-gff', dest='geneAnnotationFile', help='input gene annotation file', type=checkGffFormat, required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)
#optional:
parser.add_argument('-l', dest='lengthTerminator', help='length of terminator, default:100', type=checkInt, nargs='?', default=100)

args = parser.parse_args()

geneAnnotationFile = args.geneAnnotationFile
predictedTerminators = args.predictedTerminators
predictedNegatives = args.predictedNegatives
outpath = args.outpath
lengthTerminator = args.lengthTerminator

l1 = lengthTerminator * 0.16666666666666664
l2 = lengthTerminator - l1

l1 =  int(math.ceil(l1))
l2 = int(math.ceil(l2))


positionPredTermBed = outpath + 'Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed'
positionPredTermBedLong = outpath + 'Distance_predictedTerminators_NO_knownTerminators_NO_genes_long.bed'


positionPredNegBed = outpath + 'Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed'
positionPredNegBedLong = outpath + 'Distance_predictedNegatives_NO_knownTerminators_NO_genes_long.bed'


positionPlot = outpath + 'DistanceToGenes.png'

#######################################################################
#######################################################################

predTermData = []
predNegData = []

stopCoordsPositiveStrand = []
startCoordsNegativeStrand = []

lengthGenome = 0
organism = ''
chrom = ''
chrom2 = ''
plasmid = ''

if "BS" in predictedTerminators:
	organism = 'B.subtilis'
	chrom = 'NC_000964.3/1-4215606'
	chrom2 = 'gi|255767013|ref|NC_000964.3|'
	plasmid = ''
	lengthGenome = 4215607
if "EF" in predictedTerminators:
	organism = 'E.faecalis'
	if 'chrom' in predictedTerminators:
		chrom = "NC_004668.1"
		chrom2 = chrom
		plasmid = 'chrom.'
		lengthGenome = 3218031
	if 'pl1' in predictedTerminators:
		chrom = "NC_004669.1"
		chrom2 = chrom
		plasmid = 'pTEF1'
		lengthGenome = 66320
	if 'pl2' in predictedTerminators:
		chrom = "NC_004671.1"
		chrom2 = chrom
		plasmid = 'pTEF2'
		lengthGenome = 57660
	if 'pl3' in predictedTerminators:
		chrom = "NC_004670.1"
		chrom2 = chrom
		plasmid = 'pTEF3'
		lengthGenome = 17963
if "LM" in predictedTerminators:
	organism = 'L.monocytogenes'
	chrom = "NC_003210.1"
	chrom2 = chrom
	plasmid = ''
	lengthGenome = 2944528
if 'SP' in predictedTerminators:
	organism = 'S.pneumoniae'
	chrom = "NC_003028.3"
	chrom2 = chrom
	plasmid = ''
	lengthGenome = 2160842

print '\n' + str(organism) + ' ' + str(plasmid)
#######################################################################
#######################################################################
# read files, log transformed:

with open(predictedTerminators, 'r') as predTerm, open(predictedNegatives, 'r') as predNeg,\
		open(geneAnnotationFile, 'r') as genes:


	for line in predTerm:
		x = float(line.split(':')[1].split('_')[1])
		y = float(line.split(':')[1].split('_')[2])
		strandTerm = line.split(':')[1].split('_')[3].rstrip()
		coordTerm = int(line.split(':')[1].split('_')[0])

		predTermData.append([x,y,strandTerm,coordTerm])


	for line in predNeg:
		x = float(line.split(':')[1].split('_')[1])
		y = float(line.split(':')[1].split('_')[2])
		strandNeg = line.split(':')[1].split('_')[3].rstrip()
		coordNeg = int(line.split(':')[1].split('_')[0])

		predNegData.append([x,y,strandNeg,coordNeg])


	for line in genes:
		if '#' not in line:
			if line.split()[2] == 'gene':
				startGene = int(line.split()[3])
				endGene = int(line.split()[4])
				strandGene = line.split()[6]

				if strandGene == '+':
					stopCoordsPositiveStrand.append(endGene)
				else:
					startCoordsNegativeStrand.append(startGene)

				

#######################################################################
#######################################################################
# search closest gene (not overlapping) on opposite strand to predicted Terminators

	closestTerm = []


	for TScount1 in predTermData:
		# TScount = [x,y,strandTerm,coordTerm]

		if TScount1[2] == '-': #TScount[2] = strandTerm
			distance1 = TScount1[3] - closestGene(stopCoordsPositiveStrand, TScount1[3])

			if distance1 > 0:	
				closestTerm.append([TScount1[0]] + [TScount1[1]] + [TScount1[3]] + [distance1] + [TScount1[2]]) #[x,y,coord,distance,strand]
		else:
			distance2 = closestGene(startCoordsNegativeStrand, TScount1[3]) - TScount1[3]
			
			if distance2 > 0:
				closestTerm.append([TScount1[0]] + [TScount1[1]] + [TScount1[3]] + [distance2] + [TScount1[2]])




# search closest gene (not overlapping) on opposite strand to predicted Negatives

	closestNeg = []


	for TScount3 in predNegData:
		# TScount = [x,y,strandNeg,coordNeg]

		if TScount3[2] == '-': #TScount[2] = strandNeg
			distance5 = TScount3[3] - closestGene(stopCoordsPositiveStrand, TScount3[3])

			if distance5 > 0:
					closestNeg.append([TScount3[0]] + [TScount3[1]] + [TScount3[3]] + [distance5] + [TScount3[2]]) #[x,y,coord,distance,strand]
		else:
			distance6 = closestGene(startCoordsNegativeStrand, TScount3[3]) - TScount3[3]

			if distance6 > 0:
				closestNeg.append([TScount3[0]] + [TScount3[1]] + [TScount3[3]] + [distance6] + [TScount3[2]])




	# sorted by ascending distance:
	sortedClosestTerm = sorted(closestTerm, key=lambda x: x[3])
	sortedClosestNeg = sorted(closestNeg, key=lambda x: x[3])

	numberOfTerminatorsUnder150 = 0
	numberOfNegativesUnder150 = 0



#######################################################################
#######################################################################	
	#write original coords and distance to closest gene to bed files

	with open(positionPredTermBed, 'w') as termbed, open(positionPredNegBed, 'w') as negbed,\
		open(positionPredTermBedLong, 'w') as termbedLong, open(positionPredNegBedLong, 'w') as negbedLong:

		for item1 in sortedClosestTerm: #[x,y,coord,distance,strand]
			if item1[3] <= 150:
				numberOfTerminatorsUnder150+=1

				termbed.write(chrom2 + '\t' + str(item1[2]) + '\t' + str(item1[2]+1) + '\t' + str(item1[3]) + ' predicted Terminators' + '\t' + str(item1[4]) + '\n')

				if item1[2]-l2 > 0 and item1[2]+l2 <= lengthGenome: 
					if item1[4] == '-':
						termbedLong.write(chrom +'\t' + str(item1[2]-l2) + '\t' + str(item1[2]+l1) + '\t' \
										+ str(item1[2]) \
										 + '_'+ str(item1[3]) +'_'+ str(item1[4])+'\n')

					if item1[4] == '+':
						termbedLong.write(chrom + '\t' + str(item1[2]-l1) + '\t' + str(item1[2]+l2) + '\t' \
										+ str(item1[2]) \
										+'_'+ str(item1[3])+'_'+ str(item1[4])+ '\n')

					if item1[4] == 'nostrand':
						termbedLong.write(chrom + '\t' + str(item1[2]-l1) + '\t' + str(item1[2]+l2) + '\t' \
										+ str(item1[2]) \
										+'_'+ str(item1[3])+'_'+ '+' + '\n')
						termbedLong.write(chrom + '\t' + str(item1[2]-l2) + '\t' + str(item1[2]+l1) + '\t' \
										+ str(item1[2]) \
										+'_'+ str(item1[3])+'_'+ '-' + '\n')



		for item2 in sortedClosestNeg:
			if item2[3] <= 150:
				numberOfNegativesUnder150+=1

				negbed.write(chrom2 + '\t' + str(item2[2]) + '\t' + str(item2[2]+1) + '\t' + str(item2[3]) + ' predicted Negatives' + '\t' + str(item2[4]) +'\n')


				if item2[2]-l2 > 0 and item2[2]+l2 <= lengthGenome: 
					if item2[4] == '-':
						negbedLong.write(chrom +'\t' + str(item2[2]-l2) + '\t' + str(item2[2]+l1) + '\t' \
										+ chrom + ':'+ str(item2[2]) \
										 + '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ str(item2[4])+'\n')

					if item2[4] == '+':
						negbedLong.write(chrom + '\t' + str(item2[2]-l1) + '\t' + str(item2[2]+l2) + '\t' \
										+ chrom + ':'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ str(item2[4])+ '\n')

					if item2[4] == 'nostrand':
						negbedLong.write(chrom + '\t' + str(item2[2]-l1) + '\t' + str(item2[2]+l2) + '\t' \
										+ chrom + ':'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ '+' + '\n')
						negbedLong.write(chrom + '\t' + str(item2[2]-l2) + '\t' + str(item2[2]+l1) + '\t' \
										+ chrom + ':'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ '-' + '\n')




	print "predicted Terminators not overlapping genes, not overlapping known terminators, up to 150 nucs downstream of gene: " + str(numberOfTerminatorsUnder150)
	print "predicted Negatives not overlapping genes, not overlapping known terminators, up to 150 nucs downstream of gene: " + str(numberOfNegativesUnder150)

	negbed.close()
	negbedLong.close()
	termbed.close()
	termbedLong.close()

#######################################################################
#######################################################################

	if sortedClosestTerm:
		npClosestTerm = np.array(sortedClosestTerm)
		xvalsPredTerm = npClosestTerm[:,0]
		xvalsPredTerm = xvalsPredTerm.astype(np.float)

		yvalsPredTerm = npClosestTerm[:,1]
		yvalsPredTerm = yvalsPredTerm.astype(np.float)

		distancePredTerm = npClosestTerm[:,3]

	npClosestNeg = np.array(sortedClosestNeg)

	xvalsPredNeg = npClosestNeg[:,0]
	xvalsPredNeg = xvalsPredNeg.astype(np.float)

	yvalsPredNeg = npClosestNeg[:,1]
	yvalsPredNeg = yvalsPredNeg.astype(np.float)

	distancePredNeg = npClosestNeg[:,3]


#######################################################################
#######################################################################
	#function for decision boundary: 

	xII = np.linspace(0, 24, 12.5)
	# m = slope, b = intercept 
	m = 0.8
	b = 0.8
	# hb = horizontal intercept, hstop = x value at which horizontal line stops 
	hb = 2.5
	hstop = 2

	yII = b+m*xII
	yII[:hstop] = hb
	

#######################################################################
#######################################################################

#  predicted terminators with distance to closest gene (all Points above boundary not overlapping genes) 
#  predicted negatives with distance to closest gene (all Points below boundary)

	plt.figure(figsize=(14.5,6),dpi=120)
	plt.suptitle(str(organism) +" "+ str(plasmid), fontsize=17)

	plt.subplot(121)
	
	oldCM = plt.get_cmap('hsv',400)
	newcolors = oldCM(np.linspace(0,1,400))
	newcolor = np.array([0.143343, 0.522773, 0.556295, 1.])
	newcolors[150:, :] = newcolor
	cm = ListedColormap(newcolors)

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1], color='orange', alpha=0.5)

	if sortedClosestTerm:

		data = plt.scatter(xvalsPredTerm, yvalsPredTerm, c=distancePredTerm, s=25, label='predicted terminators', marker='P', cmap=cm, vmin = 0, vmax=400)

		cbar = plt.colorbar(data)
		cbar.ax.set_ylabel('distance to closest gene', fontsize=14)

		plt.ylim(0,12.5)
		plt.xlim(0,12.5)

		plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)', fontsize=14)
		plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)', fontsize=14)
		# plt.title('Predicted Terminators' + '\n' + '(not overl. genes, not overl. known terminators)', fontsize=14)


		plt.legend(prop={'size': 14})
		ax = plt.gca()
		legend = ax.get_legend()
		legend.legendHandles[0].set_color(plt.cm.Greys(.8))


	plt.subplot(122)


	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1], color='orange', alpha=0.5)

	data = plt.scatter(xvalsPredNeg, yvalsPredNeg, c=distancePredNeg, s=25, label='predicted negatives', marker='P', cmap=cm, vmin = 0, vmax=400)

	cbar = plt.colorbar(data)
	cbar.ax.set_ylabel('distance to closest gene', fontsize=14)

	plt.ylim(0,12.5)
	plt.xlim(0,12.5)

	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)', fontsize=14)
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)', fontsize=14)
	# plt.title('Predicted Negatives' + '\n' + '(not overl. genes, not overl. known terminators)', fontsize=14)


	plt.legend(prop={'size': 14})
	ax = plt.gca()
	legend = ax.get_legend()
	legend.legendHandles[0].set_color(plt.cm.Greys(.8))


 	plt.savefig(positionPlot,dpi=300)