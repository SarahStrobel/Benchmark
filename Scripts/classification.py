import numpy as np
import math
import matplotlib.pyplot as plt
import argparse
import os.path
from tabulate import tabulate

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
# recursive binary search for looking for overlaps with genes / known terminators etc.
# Returns index of x in arr if present, else -1 

def binarySearch (arr, l, r, x): 

    if r >= l:  
        mid = l + (r - l)/2

        if arr[mid] == x: 
            return mid 

        elif arr[mid] > x: 
            return binarySearch(arr, l, mid-1, x) 

        else: 
            return binarySearch(arr, mid+1, r, x) 
  
    else: 
        return -1

#######################################################################
#######################################################################
#functions for sens, ppv, mcc

def calcTPR(tp, fn):
    denom = float(tp + fn)
    sens = 0.00
    if (denom > 0):
        sens = tp / denom
        if sens > 1.0:
            sens = 1.0
    
    return sens

def calcFPR(fp, tn):
    denom = float(fp + tn)
    fpr = 0.00
    if (denom > 0):
        fpr = fp / denom
        if fpr > 1.0:
            fpr = 1.0
    
    return fpr

def calcPpv(tp, fp):
    denom = float(tp + fp)
    ppv = 0.00
    if (denom > 0):
        ppv = tp / denom
        if ppv > 1.0:
            ppv = 1.0
    return ppv

def calcMcc(tp, tn, fp, fn):
    denom = float(math.sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn)))
    mcc = 0.00
    if (denom > 0):
        mcc = (tp*tn - fp*fn)/denom
        if mcc > 1.0:
            mcc = 1.0
        if mcc < -1.0:
            mcc + -1.0
    return mcc


#######################################################################
#######################################################################
# calculates true positives (tn), false positives (fn), true negatives (tn), false negatives (fn), 
# true positives rate (TPR, sensitivity, recall), false positves rate (FPR) and positive predictive value (PPV, precision)
def performance(classesTest, classesPredicted):
	tp = 0
	fp = 0
	tn = 0
	fn = 0

	for i, j in zip(classesTest, classesPredicted):
		if i == 1 and j == 1:
			tp += 1	
		if i == 0 and j == 0:
			tn += 1
		if i == 0 and j == 1:
			fp += 1
		if i == 1 and j == 0:
			fn += 1

	tpr = calcTPR(tp, fn)
	fpr = calcFPR(fp, tn)
	ppv = calcPpv(tp, fp)

	return tp, fp, tn, fn, tpr, fpr, ppv



#######################################################################
#######################################################################

parser = argparse.ArgumentParser(description= 'Set boundary for classification into predicted Terminators and predicted Negatives' + '\n'
								'Usage:' + '\t' + 'classifcation.py <options> -pos -neg -all -gff -term -o -l')

#required files:
parser.add_argument('-pos', dest='positivesFile', help='input of Term-Seq counts overlapping known Terminators/counts with RNIE scores > 20.0', required=True)
parser.add_argument('-neg', dest='negativesFile', help='input of Term-Seq counts overlapping Genes', required=True)
parser.add_argument('-all', dest='allPointsFile', help='input of max. Term-Seq vs. avg. RNA-Seq counts', required=True)
parser.add_argument('-gff', dest='geneAnnotationFile', help='input gene annotation file', type=checkGffFormat, required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)
parser.add_argument('-l', dest='lengthTerminator', help='length of terminator, default:120', type=checkInt, nargs='?', default=120)

#optional (only B.subtilis)
parser.add_argument('-term', dest='knownTerminators', help='input known Terminators', type=checkBedFormat)

args = parser.parse_args()

positivesFile = args.positivesFile
negativesFile = args.negativesFile
allPointsFile = args.allPointsFile
geneAnnotationFile = args.geneAnnotationFile
outpath = args.outpath
lengthTerminator = args.lengthTerminator

l1 = lengthTerminator * 0.16666666666666664
l2 = lengthTerminator - l1

l1 =  int(math.ceil(l1))
l2 = int(math.ceil(l2))




organism = ''
organismBrev = ''
chrom = ''
chrom2 = ''
plasmid = ''
lengthGenome = 0

if 'BS' in positivesFile:
	knownTerminators = args.knownTerminators
	organism = 'B.subtilis'
	organismBrev = 'BS'
	chrom = 'NC_000964.3/1-4215606'
	chrom2 = 'gi|255767013|ref|NC_000964.3|'
	lengthGenome = 4215606
if 'EF' in positivesFile:
	organism = 'E.faecalis'
	if 'chrom' in positivesFile:
		chrom = "NC_004668.1"
		chrom2 = chrom
		organismBrev = 'EF_chrom'
		plasmid = ' chromosome'
		lengthGenome = 3218031
	if 'pl1' in positivesFile:
		chrom = "NC_004669.1"
		chrom2 = chrom
		organismBrev = 'EF_pl1'
		plasmid = ' plasmid 1'
		lengthGenome = 66320
	if 'pl2' in positivesFile:
		chrom = "NC_004671.1"
		chrom2 = chrom
		organismBrev = 'EF_pl2'
		plasmid = ' plasmid 2'
		lengthGenome = 57660
	if 'pl3' in positivesFile:
		chrom = "NC_004670.1"
		chrom2 = chrom
		organismBrev = 'EF_pl3'
		plasmid = ' plasmid 3'
		lengthGenome = 17963
if 'LM' in positivesFile:
	organism = 'L.monocytogenes'
	organismBrev = 'LM'
	chrom = "NC_003210.1"
	chrom2 = chrom
	lengthGenome = 2160842
# if 'SP' in positivesFile:
# 	organism = 'S.pneumoniae'
# 	organismBrev = 'SP'
# 	chrom = "NC_003028.3"
# 	chrom2 = chrom


print organism + plasmid
print chrom

# outfiles for plots
outfileG1 = outpath + 'gradStudentsAlgo_traintest_' + organismBrev
outfileG2 = outpath + 'gradStudentsAlgo_allPoints_' + organismBrev


# outfiles for predicted terminators bed (e.g. for IGV regions manager)/rnie bed (lengthTerminator)
allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesBed = outpath + organismBrev + '_predictedTerminators_NO_knownTerminators_NO_genes.bed'
allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesLong = outpath + organismBrev + '_predictedTerminators_NO_knownTerminators_NO_genes_long.bed'

# outfiles for predicted negatives bed/rnie bed
allPointsPredictedNegativesNOknownTerminatorsNOgenesBed = outpath + organismBrev + '_predictedNegatives_NO_knownTerminators_NO_genes.bed'
allPointsPredictedNegativesNOknownTerminatorsNOgenesLong = outpath + organismBrev + '_predictedNegatives_NO_knownTerminators_NO_genes_long.bed'


# outfile for false Positives
fpBed = outpath + organismBrev + '_falsePositives.bed'




#######################################################################
#######################################################################

testPositivesData = []
testNegativesData = []
testData = [] # testPositivesData + testNegativesData
testPositivesCoords = []
testNegativesCoords = []

testAllPointsData = []
testAllPointsCoords = []
testAllPointsStrands = []

geneCoords = []
terminatorCoords = []


#######################################################################
#######################################################################
# read files, log transformed:
if organism == 'B.subtilis':
	with open(knownTerminators, 'r') as known:
		for line in known:
			startTerminator = int(line.split()[1])+1 
			endTerminator = int(line.split()[2])+1			
			for i in range (startTerminator, endTerminator):
				terminatorCoords.append(i)

with open(positivesFile, 'r') as pf, open(negativesFile, 'r') as nf, open(allPointsFile, 'r') as apf,\
		open(geneAnnotationFile, 'r') as gene:

	for line in pf:
		x1 = float(line.split()[2])
		y1 = float(line.split()[1])
		coord1 = int(line.split()[0])
		strand1 =line.split()[3]
		testPositivesData.append([np.log(x1+1), np.log(y1+1)])
		testPositivesCoords.append(coord1)

	for line in nf:
		x2 = float(line.split()[2])
		y2 = float(line.split()[1])
		coord2 = int(line.split()[0])
		strand2 =line.split()[3]
		testNegativesData.append([np.log(x2+1), np.log(y2+1)])
		testNegativesCoords.append(coord2)

	for line in apf:
		x3 = float(line.split()[2])
		y3 = float(line.split()[1])
		coord3 = int(line.split()[0])
		strand3 =line.split()[3]
		testAllPointsData.append([np.log(x3+1), np.log(y3+1)])
		testAllPointsCoords.append(coord3)
		testAllPointsStrands.append(strand3)


	for line in gene:
		if '#' not in line:
			if line.split()[2] == 'gene':
				start = int(line.split()[3])
				end = int(line.split()[4])
				for i in range (start, end+1): #half open!
					geneCoords.append(i)




	# overlapsBoth = set(testPositivesCoords) & set(testNegativesCoords)

	# print 'known terminators overlapping genes: ' + str(overlapsBoth)

#######################################################################
#######################################################################

	npTestPositivesData = np.array(testPositivesData)
	npTestNegativesData = np.array(testNegativesData)

	testData = testPositivesData + testNegativesData

	classesTest = []
	classesTest.extend([1]*len(testPositivesData))
	classesTest.extend([0]*len(testNegativesData))


	#known classes
	#all negatives, X: datapoints, y: classes
	Xtest = np.array(testData)
	yTest = np.array(classesTest)

	#unkown classes
	XtestAllPoints = np.array(testAllPointsData)


#######################################################################
#######################################################################
	# function for own decision boundary: 

	xII = np.linspace(0, 24, 12.5)
	# m = slope, b = intercept of diagonal line
	m = 0.8
	b = 0.8
	# hb = horizontal intercept, hstop = x value at which horizontal line stops 
	hb = 2.5
	hstop = 2

	yII = b+m*xII
	yII[:hstop] = hb
	

	# test data (only overlapping genes or overlapping terminators)
	xvals = Xtest[:,0]
	yvals = Xtest[:,1]

	# test data (all points)
	xvalsAllPoints = XtestAllPoints[:,0]
	yvalsAllPoints = XtestAllPoints[:,1]


	yIII = b+m*xvals
	yIV = b+m*xvalsAllPoints


	# masks for classifying test data into positives/negatives in plots
	mask = ((yvals > yIII) & (yvals > hb))
	maskII = ((yvalsAllPoints > yIV) & (yvalsAllPoints > hb))



#######################################################################
#######################################################################
	# all values > own decision boundary = positives, all values < own decision boundary = negatives
	# get all coords to points that are predicted terminators/negatives

	# test data
	predictedTerminatorsDataTest = [] # above boundary
	predictedNegativesDataTest = [] # below boundary
	classesPredictedTest = []

	# coords for known classes test, training, validation data
	testCoords = testPositivesCoords + testNegativesCoords

	# coords for predicted classes test data
	predictedTerminatorsCoordsTest = []
	predictedNegativesCoordsTest = []	

	for x, y, coord in zip(xvals, yvals, testCoords):
		if ((y > b+m*x) & (y > hb)):
			predictedTerminatorsDataTest.append([x,y])
			classesPredictedTest.append(1)
			predictedTerminatorsCoordsTest.append(coord)
		else:
			predictedNegativesDataTest.append([x,y])
			classesPredictedTest.append(0)
			predictedNegativesCoordsTest.append(coord)


#######################################################################
#######################################################################

	# tp, fp, tn, fn, tpr, fpr, ppv for own boundary 
	tpG, fpG, tnG, fnG, tprG, fprG, ppvG = performance(classesTest, classesPredictedTest)


#######################################################################
#######################################################################

	# test data all points
	predictedTerminatorsAllPointsDataTest = [] # above boundary
	predictedNegativesAllPointsDataTest = [] # below boundary

	for x, y, coord, strand in zip(xvalsAllPoints, yvalsAllPoints, testAllPointsCoords, testAllPointsStrands):
		if ((y > b+m*x) & (y > hb)):
			predictedTerminatorsAllPointsDataTest.append([x,y,coord, strand])

		else:
			predictedNegativesAllPointsDataTest.append([x,y,coord, strand])


#######################################################################
#######################################################################
# binary search: looking for predicted terminators/negatives that overlap or don't overlap known terminators or genes in test data


	predictedTerminatorsAllPointsNOTOverlappingTest = []
	predictedTerminatorsAllPointsNOTOverlappingGenesTest = []
	predictedTerminatorsAllPointsNOTOverlappingTerminatorsTest = []

	predictedTerminatorsAllPointsOverlappingTest = []
	predictedTerminatorsAllPointsOverlappingGenesTest = []
	predictedTerminatorsAllPointsOverlappingTerminatorsTest = []


	predictedNegativesAllPointsNOTOverlappingTest = []
	predictedNegativesAllPointsNOTOverlappingGenesTest = []
	predictedNegativesAllPointsNOTOverlappingTerminatorsTest = []

	predictedNegativesAllPointsOverlappingTest = []
	predictedNegativesAllPointsOverlappingGenesTest = []
	predictedNegativesAllPointsOverlappingTerminatorsTest = []

	
	outfilesPredictedTerminators = [allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesBed, allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesLong,
									fpBed]

	outfilesPredictedNegatives = [allPointsPredictedNegativesNOknownTerminatorsNOgenesBed, allPointsPredictedNegativesNOknownTerminatorsNOgenesLong]



	outfilePredTerm = [open(i, 'w') for i in outfilesPredictedTerminators]
	outfilePredNeg = [open(i, 'w') for i in outfilesPredictedNegatives]


	for posTest in predictedTerminatorsAllPointsDataTest:
			
		# check if coord of predicted terminator overlaps with known terminators (resultPos1) of genes (resultPos2) via binary search
		resultPos1 = binarySearch(terminatorCoords, 0, len(terminatorCoords)-1, posTest[2])
		resultPos2 = binarySearch(geneCoords, 0, len(geneCoords)-1, posTest[2])
		

		# bed files with predicted terminators/predicted negatives not overl. genes, not overl. known terminators
		if resultPos1 == -1 and resultPos2 == -1:
			predictedTerminatorsAllPointsNOTOverlappingTest.append([posTest[0],posTest[1]])

			outfilePredTerm[0].write(chrom2 + '\t' + str(posTest[2]) + '\t' + str(posTest[2]+1) \
							+ '\t' + 'predicted Terminator (not overl. genes and known terminators)' + '\n')

			# make artificial terminators of 'lengthTerminator'
			# if read on negative strand: artificial terminator coords: write 'l2' upstream and 'l1' downstream from TS coord 
			# if read on positive strand: artificial terminator coords: write 'l1' upstream and 'l2' downstream from TS coord
			# if no strand: create both  
			# bed format --> bedtools getFasta

			if posTest[2]-l2 > 0 and posTest[2]+l1 <= lengthGenome: 
				if posTest[3] == '-':
					outfilePredTerm[1].write(chrom +'\t' + str(posTest[2]-l2) + '\t' + str(posTest[2]+l1) + '\t' \
									+ chrom + ':'+ str(posTest[2]) \
									 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')

				if posTest[3] == '+':
					outfilePredTerm[1].write(chrom +'\t' + str(posTest[2]-l1) + '\t' + str(posTest[2]+l2) + '\t' \
									+ chrom + ':'+ str(posTest[2]) \
									 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')
				if posTest[3] == 'nostrand':
					outfilePredTerm[1].write(chrom + '\t' + str(posTest[2]-l1) + '\t' + str(posTest[2]+l2) + '\t' \
									+ chrom + ':'+ str(posTest[2]) \
									+ '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '+' + '\n')
					outfilePredTerm[1].write(chrom +'\t' + str(posTest[2]-l2) + '\t' + str(posTest[2]+l1) + '\t' \
									+ chrom + ':'+ str(posTest[2]) \
									 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '-' +'\n')

		else:
			predictedTerminatorsAllPointsOverlappingTest.append([posTest[0],posTest[1]])

		# bed file with false positives
		if resultPos2 != -1:
			outfilePredTerm[2].write(chrom2 + '\t' + str(posTest[2]) + '\t' + str(posTest[2]+1) + '\t' + 'false positives' +'\n')


		# predicted terminators not overl. genes
		if resultPos2 == -1:
			predictedTerminatorsAllPointsNOTOverlappingGenesTest.append([posTest[0],posTest[1]])

		# predicted terminators overl. genes
		else:
			predictedTerminatorsAllPointsOverlappingGenesTest.append([posTest[0],posTest[1]])


		# predicted terminators not overl. known terminators
		if resultPos1 == -1:
			predictedTerminatorsAllPointsNOTOverlappingTerminatorsTest.append([posTest[0],posTest[1]])

		# predicted terminators overl. known terminators
		else:
			predictedTerminatorsAllPointsOverlappingTerminatorsTest.append([posTest[0],posTest[1]])




	for negTest in predictedNegativesAllPointsDataTest:

		
		resultNeg1 = binarySearch(terminatorCoords, 0, len(terminatorCoords)-1, negTest[2])
		resultNeg2 = binarySearch(geneCoords, 0, len(geneCoords)-1, negTest[2])

		if resultNeg1 == -1 and resultNeg2 ==-1:
			predictedNegativesAllPointsNOTOverlappingTest.append([negTest[0],negTest[1]])

			outfilePredNeg[0].write(chrom2 + '\t' + str(negTest[2]) + '\t' + str(negTest[2]+1) \
								+ '\t' + 'predicted negative (not overl. genes and known terminators)' + '\n')

			if negTest[2]-l2 > 0 and negTest[2]+l1 <= lengthGenome: 
				if negTest[3] == '-':
					outfilePredNeg[1].write(chrom +'\t' + str(negTest[2]-l2) + '\t' + str(negTest[2]+l1) + '\t' \
									+ chrom + ':'+ str(negTest[2]) \
									 + '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+'\n')

				if negTest[3] == '+':
					outfilePredNeg[1].write(chrom +'\t' + str(negTest[2]-l1) + '\t' + str(negTest[2]+l2) + '\t' \
									+ chrom + ':'+ str(negTest[2]) \
									 + '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+'\n')
				if negTest[3] == 'nostrand':
					outfilePredNeg[1].write(chrom + '\t' + str(negTest[2]-l1) + '\t' + str(negTest[2]+l2) + '\t' \
									+ chrom + ':'+ str(negTest[2]) \
									+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '+' + '\n')
					outfilePredNeg[1].write(chrom +'\t' + str(negTest[2]-l2) + '\t' + str(negTest[2]+l1) + '\t' \
									+ chrom + ':'+ str(negTest[2]) \
									 + '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '-' +'\n')
		else:
			predictedNegativesAllPointsOverlappingTest.append([negTest[0],negTest[1]])

		# predicted negatives not overl. genes
		if resultNeg2 == -1:
			predictedNegativesAllPointsNOTOverlappingGenesTest.append([negTest[0],negTest[1]])

		# predicted negatives overl. genes
		else:
			predictedNegativesAllPointsOverlappingGenesTest.append([negTest[0],negTest[1]])


		# predicted negatives not overl. known terminators 
		if resultNeg1 == -1:
			predictedNegativesAllPointsNOTOverlappingTerminatorsTest.append([negTest[0],negTest[1]])

		# predicted negatives overl. known terminators
		else:
			predictedNegativesAllPointsOverlappingTerminatorsTest.append([negTest[0],negTest[1]])

	
	npPredictedTerminatorsAllPointsNOTOverlappingTest = np.array(predictedTerminatorsAllPointsNOTOverlappingTest)
	npPredictedTerminatorsAllPointsOverlappingTest = np.array(predictedTerminatorsAllPointsOverlappingTest)



	print '\npredicted Terminators: ' + str(len(predictedTerminatorsAllPointsDataTest))
	print 'predicted Terminators NOT overlapping genes, NOT overlapping known terminators: ' + str(len(npPredictedTerminatorsAllPointsNOTOverlappingTest))
	print 'predicted Terminators NOT overlapping genes: ' + str(len(predictedTerminatorsAllPointsNOTOverlappingGenesTest))
	print 'predicted Terminators NOT overlapping known terminators: ' + str(len(predictedTerminatorsAllPointsNOTOverlappingTerminatorsTest))

	print 'predicted Terminators overlapping genes: ' + str(len(predictedTerminatorsAllPointsOverlappingGenesTest))
	print 'predicted Terminators overlapping known terminators: ' + str(len(predictedTerminatorsAllPointsOverlappingTerminatorsTest))


	print '\npredicted Negatives: ' + str(len(predictedNegativesAllPointsDataTest))
	print 'predicted Negatives NOT overlapping genes, NOT overlapping known terminators: ' + str(len(predictedNegativesAllPointsNOTOverlappingTest))
	print 'predicted Negatives NOT overlapping genes: ' + str(len(predictedNegativesAllPointsNOTOverlappingGenesTest))
	print 'predicted Negatives NOT overlapping known terminators: ' + str(len(predictedNegativesAllPointsNOTOverlappingTerminatorsTest))

	print 'predicted Negatives overlapping genes: ' + str(len(predictedNegativesAllPointsOverlappingGenesTest))
	print 'predicted Negatives overlapping known terminators: ' + str(len(predictedNegativesAllPointsOverlappingTerminatorsTest))





#######################################################################
#######################################################################
	# print table for precision(PPV), recall(TPR), false positives rate (FPR), TP, FP, TN, FN 

	printList =[]
	printList.append(['own boundary',ppvG, tprG, fprG, tpG, fpG, tnG, fnG])


	print '\n'
	print tabulate(printList, headers=['algorithm','PPV','TPR', 'FPR','TP', 'FP', 'TN', 'FN'], floatfmt=[".3f", ".4f", ".4f", ".4f", ".0f", ".0f", ".0f", ".0f"])
	print '\n'



#######################################################################
#######################################################################
# plot own decision boundary:
#############################
	


	fig = plt.figure(figsize=(12,12),dpi=120)

	plt.suptitle('Own Decision boundary, PPV: ' + str(ppvG) + ', intercept (b): ' + str(b) + ', slope (m): ' + str(m))



#####################################################################################


	ax = fig.add_subplot(2,2,1)
	ax.text(-0.1, 1.1, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='orange', alpha=0.5)

	classes = yTest

	unique = list(set(classes))
	colors = ['rebeccapurple', '#EE6C00']

	labelsPredicted = ['overlapping genes', 'RNIE score > 20.0']
	if organism == 'B.subtilis':
		labelsPredicted = ['overlapping genes', 'overlapping terminators']

	for i, u in enumerate(unique):
	    xi = [Xtest[:,0][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]
	    yi = [Xtest[:,1][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, label=str(labelsPredicted[i]), marker='x')

	plt.ylim(0,12.5)
	plt.xlim(0,12.5)
	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)')
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)')
	plt.title('known classes')
	plt.legend()	

#####################################################################################


	ax = fig.add_subplot(2,2,2)
	ax.text(-0.1, 1.1, 'B', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='orange', alpha=0.5)


	plt.scatter(xvals[~mask], yvals[~mask], c='#8E5BE9', s=12, label='predicted negatives', marker='x')
	plt.scatter(xvals[mask], yvals[mask], c='#FFA913', s=12, label='predicted positives', marker='x')

	plt.ylim(0,12.5)
	plt.xlim(0,12.5)
	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)')
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)')
	plt.title('predicted classes')
	plt.legend()




#####################################################################################


	ax = fig.add_subplot(2,2,3)
	ax.text(-0.1, 1.1, 'C', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='orange', alpha=0.5)

	masknew = classes==0
	classesNew = classes[masknew]
	XtestNew = Xtest[masknew]
	unique = list(set(classesNew))
	colors = ['rebeccapurple']

	labelsPredicted = ['overlapping genes']
	if organism == 'B.subtilis':
		labelsPredicted = ['overlapping genes']

	for i, u in enumerate(unique):
	    xi = [XtestNew[:,0][j] for j  in range(len(XtestNew[:,0])) if classesNew[j] == u]
	    yi = [XtestNew[:,1][j] for j  in range(len(XtestNew[:,0])) if classesNew[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, label=str(labelsPredicted[i]), marker='x')

	plt.ylim(0,12.5)
	plt.xlim(0,12.5)
	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)')
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)')
	plt.title('known classes')
	plt.legend()

#####################################################################################


	ax = fig.add_subplot(2,2,4)
	ax.text(-0.1, 1.1, 'D', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='orange', alpha=0.5)

	masknew = classes==1
	classesNew = classes[masknew]
	XtestNew = Xtest[masknew]

	unique = list(set(classesNew))
	colors = ['#EE6C00']

	labelsPredicted = ['RNIE score > 20.0']
	if organism == 'B.subtilis':
		labelsPredicted = ['overlapping terminators']

	for i, u in enumerate(unique):
	    xi = [XtestNew[:,0][j] for j  in range(len(XtestNew[:,0])) if classesNew[j] == u]
	    yi = [XtestNew[:,1][j] for j  in range(len(XtestNew[:,0])) if classesNew[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, label=str(labelsPredicted[i]), marker='x')

	plt.ylim(0,12.5)
	plt.xlim(0,12.5)
	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)')
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)')
	plt.title('known classes')
	plt.legend()

	plt.savefig(outfileG1,dpi=300)



####################################################################################
####################################################################################

	fig2 = plt.figure(figsize=(12,12),dpi=120)
	plt.suptitle('Own Decision boundary, PPV: ' + str(ppvG) + ', intercept (b): ' + str(b) + ', slope (m): ' + str(m))

	unique = list(set(classes))
	colors = ['rebeccapurple', '#EE6C00']

	labelsPredicted = ['overlapping genes', 'RNIE score > 20.0']
	if organism == 'B.subtilis':
		labelsPredicted = ['overlapping genes', 'overlapping terminators']

	#plot known classes, predicted classes
	ax = fig2.add_subplot(2,2,1)
	ax.text(-0.1, 1.1, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='orange', alpha=0.5)


	plt.scatter(xvalsAllPoints[~maskII], yvalsAllPoints[~maskII], c='#8E5BE9', s=12, label='predicted negatives', marker='+')
	plt.scatter(xvalsAllPoints[maskII], yvalsAllPoints[maskII], c='#FFA913', s=12, label='predicted positives', marker='+')


	for i, u in enumerate(unique):
	    xi = [Xtest[:,0][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]
	    yi = [Xtest[:,1][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, label=str(labelsPredicted[i]), marker='x')


	plt.ylim(0,12.5)
	plt.xlim(0,12.5)

	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)')
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)')


	plt.legend(bbox_to_anchor=(1.04,0.5), loc='center left')

################################################################################
	#plot known classes, predicted classes, overlapping known classes (red crosses)

	ax = fig2.add_subplot(2,2,3)
	ax.text(-0.1, 1.1, 'B', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='orange', alpha=0.5)


	plt.scatter(xvalsAllPoints[~maskII], yvalsAllPoints[~maskII], c='#8E5BE9', s=12, label='predicted negatives', marker='+')
	plt.scatter(xvalsAllPoints[maskII], yvalsAllPoints[maskII], c='#FFA913', s=12, label='predicted positives', marker='+')



	if not npPredictedTerminatorsAllPointsOverlappingTest.size == 0:
		plt.scatter(npPredictedTerminatorsAllPointsOverlappingTest[:,0],npPredictedTerminatorsAllPointsOverlappingTest[:,1], c='red',s=11, label='overlapping', marker ='+')
	else:
		plt.scatter(0,0, c='red',s=11, label='overlapping', marker ='+' )

	plt.ylim(0,12.5)
	plt.xlim(0,12.5)

	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)')
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)')


	plt.legend()

##############################################################################
	#plot known classes, predicted classes, NOT overlapping known classes (red crosses)

	ax = fig2.add_subplot(2,2,4)
	ax.text(-0.1, 1.1, 'C', transform=ax.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='mediumpurple', alpha=0.5)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='orange', alpha=0.5)


	plt.scatter(xvalsAllPoints[~maskII], yvalsAllPoints[~maskII], c='#8E5BE9', s=12, label='predicted negatives', marker='+')
	plt.scatter(xvalsAllPoints[maskII], yvalsAllPoints[maskII], c='#FFA913', s=12, label='predicted positives', marker='+')

	if not npPredictedTerminatorsAllPointsNOTOverlappingTest.size == 0:
		plt.scatter(npPredictedTerminatorsAllPointsNOTOverlappingTest[:,0],npPredictedTerminatorsAllPointsNOTOverlappingTest[:,1], c='red',s=11, label='not overlapping', marker ='+')
	else:
		plt.scatter(0,0, c='red',s=11, label='not overlapping', marker ='+' )

	plt.ylim(0,12.5)
	plt.xlim(0,12.5)

	plt.xlabel('$log_{e}$' + '(Avg. RNA-Seq count)')
	plt.ylabel('$log_{e}$' + '(Max. Term-Seq count)')


	plt.legend()

	plt.savefig(outfileG2,dpi=300)


