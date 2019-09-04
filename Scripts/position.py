# takes output from ratioTermRNASeq_replicates_500k.py and makes scatterplots of training data +
# takes output from SVM_gradStudAlgo_combined.py (bed files) and adds scatterplots of  
# predicted Terminators and predicted Negatives + colors according to their distances to annotated gene


import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect_left
from matplotlib.colors import ListedColormap
import argparse

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
								'Usage:' + '\t' + 'position.py <options> -pos -neg -all -gene -term -o')

#required files:
parser.add_argument('-pos', dest='predictedTerminators', help='input predicted Terminators', required=True)
parser.add_argument('-neg', dest='predictedNegatives', help='input predicted Negatives', required=True)
parser.add_argument('-gene', dest='geneAnnotationFile', help='input gene annotation file', required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)

args = parser.parse_args()

geneAnnotationFile = args.geneAnnotationFile
predictedTerminators = args.predictedTerminators
predictedNegatives = args.predictedNegatives
outpath = args.outpath


positionPredTermBed = outpath + 'Distance_predictedTerminators_NO_genes.bed'
positionPredTermBed120 = outpath + 'Distance_120_predictedTerminators_NO_genes.bed'
positionPredTermBed60 = outpath + 'Distance_60_predictedTerminators_NO_genes.bed'

positionPredNegBed = outpath + 'Distance_predictedNegatives_NO_genes.bed'
positionPredNegBed120 = outpath + 'Distance_120_predictedNegatives_NO_genes.bed'
positionPredNegBed60 = outpath + 'Distance_60_predictedNegatives_NO_genes.bed'

positionPlot = outpath + 'DistanceToGenes.png'

#######################################################################
#######################################################################

predTermData = []
predNegData = []

stopCoordsPositiveStrand = []
startCoordsNegativeStrand = []


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
# search closest gene (not overlapping) on opposite strand to predictedTerminators

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




# search closest gene (not overlapping) on opposite strand to predictedNegatives

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





#######################################################################
#######################################################################	
	#write original coords and distance to closest gene to bed files

	with open(positionPredTermBed, 'w') as termbed, open(positionPredNegBed, 'w') as negbed,\
		open(positionPredTermBed120, 'w') as termbed120, open(positionPredNegBed120, 'w') as negbed120,\
		open(positionPredTermBed60, 'w') as termbed60, open(positionPredNegBed60, 'w') as negbed60:

		for item1 in sortedClosestTerm: #[x,y,coord,distance,strand]
			if item1[3] <= 150:
				termbed.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(item1[2]) + '\t' + str(item1[2]+1) + '\t' + str(item1[3]) + ' predicted Terminators' + '\t' + str(item1[4]) + '\n')

				if item1[2]-100 > 0 and item1[2]+100 <= 4215606: 
					if item1[4] == '-':
						termbed120.write('NC_000964.3/1-4215606' +'\t' + str(item1[2]-100) + '\t' + str(item1[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										 + '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3]) +'_'+ str(item1[4])+'\n')

					if item1[4] == '+':
						termbed120.write('NC_000964.3/1-4215606' + '\t' + str(item1[2]-20) + '\t' + str(item1[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										+ '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3])+'_'+ str(item1[4])+ '\n')

					if item1[4] == 'nostrand':
						termbed120.write('NC_000964.3/1-4215606' + '\t' + str(item1[2]-20) + '\t' + str(item1[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										+ '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3])+'_'+ '+' + '\n')
						termbed120.write('NC_000964.3/1-4215606' + '\t' + str(item1[2]-100) + '\t' + str(item1[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										+ '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3])+'_'+ '-' + '\n')

				if item1[2]-50 > 0 and item1[2]+50 <= 4215606: 
					if item1[4] == '-':
						termbed60.write('NC_000964.3/1-4215606' +'\t' + str(item1[2]-50) + '\t' + str(item1[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										 + '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3])+'_'+ str(item1[4])+'\n')

					if item1[4] == '+':
						termbed60.write('NC_000964.3/1-4215606' + '\t' + str(item1[2]-10) + '\t' + str(item1[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										+ '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3])+'_'+ str(item1[4])+ '\n')

					if item1[4] == 'nostrand':
						termbed60.write('NC_000964.3/1-4215606' + '\t' + str(item1[2]-10) + '\t' + str(item1[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										+ '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3])+'_'+ '+' + '\n')
						termbed60.write('NC_000964.3/1-4215606' + '\t' + str(item1[2]-50) + '\t' + str(item1[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item1[2]) \
										+ '_' + str(item1[0]) + '_' + str(item1[1]) +'_'+ str(item1[3])+'_'+ '-' + '\n')





		for item2 in sortedClosestNeg:
			if item2[3] <= 150:
				negbed.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(item2[2]) + '\t' + str(item2[2]+1) + '\t' + str(item2[3]) + ' predicted Negatives' + '\t' + str(item2[4]) +'\n')


				if item2[2]-100 > 0 and item2[2]+100 <= 4215606: 
					if item2[4] == '-':
						negbed120.write('NC_000964.3/1-4215606' +'\t' + str(item2[2]-100) + '\t' + str(item2[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										 + '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ str(item2[4])+'\n')

					if item2[4] == '+':
						negbed120.write('NC_000964.3/1-4215606' + '\t' + str(item2[2]-20) + '\t' + str(item2[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ str(item2[4])+ '\n')

					if item2[4] == 'nostrand':
						negbed120.write('NC_000964.3/1-4215606' + '\t' + str(item2[2]-20) + '\t' + str(item2[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ '+' + '\n')
						negbed120.write('NC_000964.3/1-4215606' + '\t' + str(item2[2]-100) + '\t' + str(item2[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ '-' + '\n')

				if item2[2]-50 > 0 and item2[2]+50 <= 4215606: 
					if item2[4] == '-':
						negbed60.write('NC_000964.3/1-4215606' +'\t' + str(item2[2]-50) + '\t' + str(item2[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										 + '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ str(item2[4])+'\n')

					if item2[4] == '+':
						negbed60.write('NC_000964.3/1-4215606' + '\t' + str(item2[2]-10) + '\t' + str(item2[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ str(item2[4])+ '\n')

					if item2[4] == 'nostrand':
						negbed60.write('NC_000964.3/1-4215606' + '\t' + str(item2[2]-10) + '\t' + str(item2[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ '+' + '\n')
						negbed60.write('NC_000964.3/1-4215606' + '\t' + str(item2[2]-50) + '\t' + str(item2[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(item2[2]) \
										+ '_' + str(item2[0]) + '_' + str(item2[1]) +'_'+ str(item2[3])+'_'+ '-' + '\n')


	negbed.close()
	negbed60.close()
	negbed120.close()
	termbed.close()
	termbed60.close()
	termbed120.close()

#######################################################################
#######################################################################

	npClosestTerm = np.array(sortedClosestTerm)
	npClosestNeg = np.array(sortedClosestNeg)


	xvalsPredTerm = npClosestTerm[:,0]
	xvalsPredTerm = xvalsPredTerm.astype(np.float)
	yvalsPredTerm = npClosestTerm[:,1]
	yvalsPredTerm = yvalsPredTerm.astype(np.float)
	distancePredTerm = npClosestTerm[:,3]

	xvalsPredNeg = npClosestNeg[:,0]
	xvalsPredNeg = xvalsPredNeg.astype(np.float)
	yvalsPredNeg = npClosestNeg[:,1]
	yvalsPredNeg = yvalsPredNeg.astype(np.float)
	distancePredNeg = npClosestNeg[:,3]


#######################################################################
#######################################################################
	#function for own decision boundary: 

	xII = np.linspace(0, 24, 12.5)
	# m = slope, b = intercept 
	m = 0.7
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
	plt.suptitle('Distance to closest gene')

	plt.subplot(121)

	oldCM = plt.get_cmap('hsv',400)
	newcolors = oldCM(np.linspace(0,1,400))
	newcolor = np.array([0.143343, 0.522773, 0.556295, 1.])
	newcolors[150:, :] = newcolor
	cm = ListedColormap(newcolors)


	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='#E87CC6', alpha=0.8)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='#7A9441', alpha=0.5)

	data = plt.scatter(xvalsPredTerm, yvalsPredTerm, c=distancePredTerm, s=25, label='predicted terminators - distance to gene', marker='P', cmap=cm, vmin = 0, vmax=400)

	cbar = plt.colorbar(data)
	cbar.ax.set_ylabel('distance to closest gene')

	plt.ylim(0,10.5)
	plt.xlim(0,12.5)

	plt.ylabel('Max. Term-Seq count')
	plt.xlabel('Avg. RNA-Seq count')
	plt.title('Predicted Terminators' + '\n' + '(not overl. genes)')


	plt.legend()
	ax = plt.gca()
	legend = ax.get_legend()
	legend.legendHandles[0].set_color(plt.cm.Greys(.8))


	plt.subplot(122)


	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='#E87CC6', alpha=0.8)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='#7A9441', alpha=0.5)

	data = plt.scatter(xvalsPredNeg, yvalsPredNeg, c=distancePredNeg, s=25, label='predicted negatives - distance to gene', marker='P', cmap=cm, vmin = 0, vmax=400)

	cbar = plt.colorbar(data)
	cbar.ax.set_ylabel('distance to closest gene')

	plt.ylim(0,10.5)
	plt.xlim(0,12.5)

	plt.ylabel('Max. Term-Seq count')
	plt.xlabel('Avg. RNA-Seq count')
	plt.title('Predicted Negatives' + '\n' + '(not overl. genes)')


	plt.legend()
	ax = plt.gca()
	legend = ax.get_legend()
	legend.legendHandles[0].set_color(plt.cm.Greys(.8))


 	plt.savefig(positionPlot,dpi=300)