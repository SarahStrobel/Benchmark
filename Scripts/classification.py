import numpy as np
import matplotlib.pyplot as plt
import argparse
from tabulate import tabulate




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
								'Usage:' + '\t' + 'classifcation.py <options> -pos -neg -all -gene -term -o')

#required files:
parser.add_argument('-pos', dest='positivesFile', help='input of Term-Seq counts overlapping known Terminators', required=True)
parser.add_argument('-neg', dest='negativesFile', help='input of Term-Seq counts overlapping Genes', required=True)
parser.add_argument('-all', dest='allPointsFile', help='input of max. Term-Seq vs. avg. RNA-Seq counts', required=True)
parser.add_argument('-gene', dest='geneAnnotationFile', help='input gene annotation file', required=True)
parser.add_argument('-term', dest='knownTerminators', help='input known Terminators', required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)

args = parser.parse_args()

positivesFile = args.positivesFile
negativesFile = args.negativesFile
allPointsFile = args.allPointsFile
geneAnnotationFile = args.geneAnnotationFile
knownTerminators = args.knownTerminators
outpath = args.outpath




outfileG1 = outpath + 'gradStudentsAlgo_traintest_'
outfileG2 = outpath + 'gradStudentsAlgo_allPoints_'



# outfiles for predicted terminators bed/rnie bed
allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesBed = outpath + 'predictedTerminators_NO_knownTerminators_NO_genes.bed'
allPointsPredictedTerminatorsNOgenesBed = outpath + 'predictedTerminators_NO_genes.bed'
allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesRNIE = outpath + '120_predictedTerminators_NO_knownTerminators_NO_genes.bed'
allPointsPredictedTerminatorsNOgenesRNIE = outpath + '120_predictedTerminators_NO_genes.bed'
allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesRNIE2 = outpath + '60_predictedTerminators_NO_knownTerminators_NO_genes.bed'
allPointsPredictedTerminatorsNOgenesRNIE2 = outpath + '60_predictedTerminators_NO_genes.bed'

# outfiles for predicted negatives bed/rnie bed
allPointsPredictedNegativesNOknownTerminatorsNOgenesBed = outpath + 'predictedNegatives_NO_knownTerminators_NO_genes.bed'
allPointsPredictedNegativesNOgenesBed = outpath + 'predictedNegatives_NO_genes.bed'
allPointsPredictedNegativesNOknownTerminatorsNOgenesRNIE = outpath + '120_predictedNegatives_NO_knownTerminators_NO_genes.bed'
allPointsPredictedNegativesNOgenesRNIE = outpath + '120_predictedNegatives_NO_genes.bed'
allPointsPredictedNegativesNOknownTerminatorsNOgenesRNIE2 = outpath + '60_predictedNegatives_NO_knownTerminators_NO_genes.bed'
allPointsPredictedNegativesNOgenesRNIE2 = outpath + '60_predictedNegatives_NO_genes.bed'

# outfile for false Positives

fpBed = outpath + 'falsePositives.bed'




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

with open(positivesFile, 'r') as pf, open(negativesFile, 'r') as nf, open(allPointsFile, 'r') as apf,\
		open(geneAnnotationFile, 'r') as gene, open(knownTerminators, 'r') as known:


	for line1 in pf:
		x1 = float(line1.split()[2])
		y1 = float(line1.split()[1])
		coord1 = int(line1.split()[0])
		strand1 =line1.split()[3]
		testPositivesData.append([np.log(x1+1), np.log(y1+1)])
		testPositivesCoords.append(coord1)


	for line2 in nf:
		x2 = float(line2.split()[2])
		y2 = float(line2.split()[1])
		coord2 = int(line2.split()[0])
		strand2 =line2.split()[3]
		testNegativesData.append([np.log(x2+1), np.log(y2+1)])
		testNegativesCoords.append(coord2)

	for line3 in apf:
		x3 = float(line3.split()[2])
		y3 = float(line3.split()[1])
		coord3 = int(line3.split()[0])
		strand3 =line3.split()[3]
		testAllPointsData.append([np.log(x3+1), np.log(y3+1)])
		testAllPointsCoords.append(coord3)
		testAllPointsStrands.append(strand3)


	for line7 in gene:
		if '#' not in line7:
			if line7.split()[2] == 'gene':
				start = int(line7.split()[3])
				end = int(line7.split()[4])
				for i in range (start, end+1): #half open!
					geneCoords.append(i)

	for line8 in known:
			startTerminator = int(line8.split()[1])+1 
			endTerminator = int(line8.split()[2])+1			
			for i in range (startTerminator, endTerminator):
				terminatorCoords.append(i)



	print 'all files read'


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
	# m = 0.8
	# b = 0.5
	m = 0.7
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
# binary search: looking for predicted positves from test data (all points) that overlap or don't overlap known terminators or genes in test data


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


	with open(allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesBed, 'w') as aPpTnOkTnOgBed, open(allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesRNIE, 'w') as aPpTnOkTnOgRNIE, \
			open(allPointsPredictedTerminatorsNOgenesBed, 'w') as aPpTnOkTbed, open(allPointsPredictedTerminatorsNOgenesRNIE, 'w')as aPpTnOkTrnie, \
			open(fpBed, 'w') as fBed,\
			open(allPointsPredictedTerminatorsNOknwonTerminatorsNOgenesRNIE2, 'w') as aPpTnOkTnOgRNIE2, \
			open(allPointsPredictedTerminatorsNOgenesRNIE2, 'w')as aPpTnOkTrnie2:


		for posTest in predictedTerminatorsAllPointsDataTest:
			
			# check if coord of predicted terminator overlaps with known terminators (resultPos1) of genes (resultPos2) via binary search
			# resultPos1 = binarySearch(testPositivesCoords, 0, len(testPositivesCoords)-1, posTest[2])
			resultPos1 = binarySearch(terminatorCoords, 0, len(terminatorCoords)-1, posTest[2])
			# resultPos2 = binarySearch(testNegativesCoords, 0, len(testNegativesCoords)-1, posTest[2])
			resultPos2 = binarySearch(geneCoords, 0, len(geneCoords)-1, posTest[2])


			# bed file and rnie bed file with predicted terminators not overl. genes, not overl. known terminators
			if resultPos1 == -1 and resultPos2 == -1:
				predictedTerminatorsAllPointsNOTOverlappingTest.append([posTest[0],posTest[1]])

				# write to bed: all points predicted Terminators Not overlapping genes or known terminators (in test data)
				aPpTnOkTnOgBed.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(posTest[2]) + '\t' + str(posTest[2]+1) \
								+ '\t' + 'predicted Terminator (not overl. genes and known terminators)' + '\n')


				# make 'artificial' terminators for RNIE
				# if read on negative strand: artificial terminator coords: write 70 upstream and 20 downstream from TS coord 
				# if read on positive strand: artificial terminator coords: write 20 upstream and 70 downstream from TS coord
				# if no strand: create both  
				# to bed format --> bedtools getFasta

				# write to bed: chrom tab coord-100/coord-20 tab coord+20/coord+100 tab name(chrom:originalTScoord_x_y_strand)

				if posTest[2]-100 > 0 and posTest[2]+100 <= 4215606: 
					if posTest[3] == '-':
						aPpTnOkTnOgRNIE.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-100) + '\t' + str(posTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')

					if posTest[3] == '+':
						aPpTnOkTnOgRNIE.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-20) + '\t' + str(posTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')
					if posTest[3] == 'nostrand':
						aPpTnOkTrnie.write('NC_000964.3/1-4215606' + '\t' + str(posTest[2]-20) + '\t' + str(posTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										+ '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '+' + '\n')
						aPpTnOkTrnie.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-100) + '\t' + str(posTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '-' +'\n')

				# write to bed: chrom tab coord-50/coord-10 tab coord+10/coord+50 tab name(chrom:originalTScoord_x_y_strand)		
				if posTest[2]-50 > 0 and posTest[2]+50 <= 4215606: 
					if posTest[3] == '-':
						aPpTnOkTnOgRNIE2.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-50) + '\t' + str(posTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')

					if posTest[3] == '+':
						aPpTnOkTnOgRNIE2.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-10) + '\t' + str(posTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')
					if posTest[3] == 'nostrand':
						aPpTnOkTrnie2.write('NC_000964.3/1-4215606' + '\t' + str(posTest[2]-10) + '\t' + str(posTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										+ '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '+' + '\n')
						aPpTnOkTrnie2.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-50) + '\t' + str(posTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '-' +'\n')


			else:
				predictedTerminatorsAllPointsOverlappingTest.append([posTest[0],posTest[1]])

			
			# bed file with false positives
			if resultPos2 != -1:
				fBed.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(posTest[2]) + '\t' + str(posTest[2]+1) + '\t' + 'false positives' +'\n')


			# bed file with predicted terminators not overl. genes alone / overlapping genes alone
			if resultPos2 == -1:
				predictedTerminatorsAllPointsNOTOverlappingGenesTest.append([posTest[0],posTest[1]])

				aPpTnOkTbed.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(posTest[2]) + '\t' + str(posTest[2]+1) \
								+ '\t' + 'predicted Terminator (not overl. genes)' + '\n')

				if posTest[2]-100 > 0 and posTest[2]+100 <= 4215606:
					if posTest[3] == '-':
						aPpTnOkTrnie.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-100) + '\t' + str(posTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')

					if posTest[3] == '+':
						aPpTnOkTrnie.write('NC_000964.3/1-4215606' + '\t' + str(posTest[2]-20) + '\t' + str(posTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										+ '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+ '\n')

					if posTest[3] == 'nostrand':
						aPpTnOkTrnie.write('NC_000964.3/1-4215606' + '\t' + str(posTest[2]-20) + '\t' + str(posTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										+ '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '+' + '\n')
						aPpTnOkTrnie.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-100) + '\t' + str(posTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '-' +'\n')


				if posTest[2]-50 > 0 and posTest[2]+50 <= 4215606:
					if posTest[3] == '-':
						aPpTnOkTrnie2.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-50) + '\t' + str(posTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+'\n')

					if posTest[3] == '+':
						aPpTnOkTrnie2.write('NC_000964.3/1-4215606' + '\t' + str(posTest[2]-10) + '\t' + str(posTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										+ '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ str(posTest[3])+ '\n')

					if posTest[3] == 'nostrand':
						aPpTnOkTrnie2.write('NC_000964.3/1-4215606' + '\t' + str(posTest[2]-10) + '\t' + str(posTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										+ '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '+' + '\n')
						aPpTnOkTrnie2.write('NC_000964.3/1-4215606' +'\t' + str(posTest[2]-50) + '\t' + str(posTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(posTest[2]) \
										 + '_' + str(posTest[0]) + '_' + str(posTest[1]) +'_'+ '-' +'\n')

			else:
				predictedTerminatorsAllPointsOverlappingGenesTest.append([posTest[0],posTest[1]])



			# predicted terminators not overl. known terminators alone / overlapping known terminators alone
			if resultPos1 == -1:
				predictedTerminatorsAllPointsNOTOverlappingTerminatorsTest.append([posTest[0],posTest[1]])

			else:
				predictedTerminatorsAllPointsOverlappingTerminatorsTest.append([posTest[0],posTest[1]])




	with open(allPointsPredictedNegativesNOknownTerminatorsNOgenesBed, 'w') as aPpNnOkTnOgBed, open(allPointsPredictedNegativesNOknownTerminatorsNOgenesRNIE, 'w') as aPpNnOkTnOgRNIE,\
			open(allPointsPredictedNegativesNOgenesBed, 'w') as aPpNnOkTbed,  open(allPointsPredictedNegativesNOgenesRNIE, 'w') as aPpNnOkTrnie, \
			open(allPointsPredictedNegativesNOknownTerminatorsNOgenesRNIE2, 'w') as aPpNnOkTnOgRNIE2,\
			open(allPointsPredictedNegativesNOgenesRNIE2, 'w') as aPpNnOkTrnie2:

		for negTest in predictedNegativesAllPointsDataTest:

			# resultNeg1 = binarySearch(testPositivesCoords, 0, len(testPositivesCoords)-1, negTest[2])
			resultNeg1 = binarySearch(terminatorCoords, 0, len(terminatorCoords)-1, negTest[2])
			# resultNeg2 = binarySearch(testNegativesCoords, 0, len(testNegativesCoords)-1, negTest[2])
			resultNeg2 = binarySearch(geneCoords, 0, len(geneCoords)-1, negTest[2])

			# bed file and rnie bed file with predicted negatives not overl. genes, not overl. known terminators
			if resultNeg1 == -1 and resultNeg2 ==-1:
				predictedNegativesAllPointsNOTOverlappingTest.append([negTest[0],negTest[1]])

				aPpNnOkTnOgBed.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(negTest[2]) + '\t' + str(negTest[2]+1) \
									+ '\t' + 'predicted negative (not overl. genes and known terminators)' + '\n')


				if negTest[2]-100 > 0 and negTest[2]+100 <= 4215606: 
					if negTest[3] == '-':
						aPpNnOkTnOgRNIE.write('NC_000964.3/1-4215606' +'\t' + str(negTest[2]-100) + '\t' + str(negTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										 + '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+'\n')

					if negTest[3] == '+':
						aPpNnOkTnOgRNIE.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-20) + '\t' + str(negTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+ '\n')

					if negTest[3] == 'nostrand':
						aPpNnOkTnOgRNIE.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-20) + '\t' + str(negTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '+' + '\n')
						aPpNnOkTnOgRNIE.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-100) + '\t' + str(negTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '-' + '\n')

				if negTest[2]-50 > 0 and negTest[2]+50 <= 4215606: 
					if negTest[3] == '-':
						aPpNnOkTnOgRNIE2.write('NC_000964.3/1-4215606' +'\t' + str(negTest[2]-50) + '\t' + str(negTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										 + '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+'\n')

					if negTest[3] == '+':
						aPpNnOkTnOgRNIE2.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-10) + '\t' + str(negTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+ '\n')

					if negTest[3] == 'nostrand':
						aPpNnOkTnOgRNIE2.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-10) + '\t' + str(negTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '+' + '\n')
						aPpNnOkTnOgRNIE2.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-50) + '\t' + str(negTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '-' + '\n')

			else:
				predictedNegativesAllPointsOverlappingTest.append([negTest[0],negTest[1]])

			# bed file with predicted negatives not overl. genes alone / overlapping genes alone
			if resultNeg2 == -1:
				predictedNegativesAllPointsNOTOverlappingGenesTest.append([negTest[0],negTest[1]])

				aPpNnOkTbed.write('gi|255767013|ref|NC_000964.3|' + '\t' + str(negTest[2]) + '\t' + str(negTest[2]+1) \
									+ '\t' + 'predicted negative (not overl. genes)' + '\n')
				if negTest[2]-100 > 0 and negTest[2]+100 <= 4215606:
					if negTest[3] == '-':
						aPpNnOkTrnie.write('NC_000964.3/1-4215606' +'\t' + str(negTest[2]-100) + '\t' + str(negTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										 + '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+'\n')

					if negTest[3] == '+':
						aPpNnOkTrnie.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-20) + '\t' + str(negTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+ '\n')

					if negTest[3] == 'nostrand':
						aPpNnOkTrnie.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-20) + '\t' + str(negTest[2]+100) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '+' + '\n')
						aPpNnOkTrnie.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-100) + '\t' + str(negTest[2]+20) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '-' + '\n')

				if negTest[2]-50 > 0 and negTest[2]+50 <= 4215606:
					if negTest[3] == '-':
						aPpNnOkTrnie2.write('NC_000964.3/1-4215606' +'\t' + str(negTest[2]-50) + '\t' + str(negTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										 + '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+'\n')

					if negTest[3] == '+':
						aPpNnOkTrnie2.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-10) + '\t' + str(negTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ str(negTest[3])+ '\n')

					if negTest[3] == 'nostrand':
						aPpNnOkTrnie2.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-10) + '\t' + str(negTest[2]+50) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '+' + '\n')
						aPpNnOkTrnie2.write('NC_000964.3/1-4215606' + '\t' + str(negTest[2]-50) + '\t' + str(negTest[2]+10) + '\t' \
										+ 'NC_000964.3/1-4215606:'+ str(negTest[2]) \
										+ '_' + str(negTest[0]) + '_' + str(negTest[1]) +'_'+ '-' + '\n')

			else:
				predictedNegativesAllPointsOverlappingGenesTest.append([negTest[0],negTest[1]])


			# predicted negatives not overl. known terminators alone / overlapping knwon terminators alone
			if resultNeg1 == -1:
				predictedNegativesAllPointsNOTOverlappingTerminatorsTest.append([negTest[0],negTest[1]])
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


	aPpNnOkTnOgBed.close()
	aPpNnOkTnOgRNIE.close()
	aPpNnOkTbed.close()
	aPpNnOkTrnie.close()	

	aPpTnOkTnOgBed.close()
	aPpTnOkTnOgRNIE.close()
	fBed.close()
	aPpTnOkTbed.close()
	aPpTnOkTrnie.close()


#######################################################################
#######################################################################
	# print table for precision(PPV), recall(TPR), false positives rate (FPR), TP, FP, TN, FN for different algos


	printList =[]
	printList.append(['own boundary',ppvG, tprG, fprG, tpG, fpG, tnG, fnG])


	print '\n'
	print tabulate(printList, headers=['algorithm','PPV','TPR', 'FPR','TP', 'FP', 'TN', 'FN'], floatfmt=[".3f", ".4f", ".4f", ".4f", ".0f", ".0f", ".0f", ".0f"])
	# print tabulate([['SVM 500', ]])
	print '\n'




#######################################################################
#######################################################################
# plot own decision boundary:
#############################





	plt.figure(figsize=(12,6),dpi=120)

	plt.suptitle('Own Decision boundary, PPV: ' + str(ppvG) + ', intercept (b): ' + str(b) + ', slope (m): ' + str(m))



#####################################################################################


	plt.subplot(121)

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='#E87CC6', alpha=0.8)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='#7A9441', alpha=0.5)

	classes = yTest


	unique = list(set(classes))
	colors = ['darkmagenta', 'darkgreen']
	labelsPredicted = ['overlapping genes', 'overlapping terminators']

	for i, u in enumerate(unique):
	    xi = [Xtest[:,0][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]
	    yi = [Xtest[:,1][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, label=str(labelsPredicted[i]), marker='x')

	plt.ylim(0,8.5)
	plt.xlim(0,12.5)
	plt.xlabel('Avg. RNA-Seq count')
	plt.ylabel('Max. Term-Seq count')
	plt.title('known classes')
	plt.legend()	

#####################################################################################


	plt.subplot(122)

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='#E87CC6', alpha=0.8)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='#7A9441', alpha=0.5)


	plt.scatter(xvals[~mask], yvals[~mask], c='darkmagenta', s=12, label='predicted negatives', marker='x')
	plt.scatter(xvals[mask], yvals[mask], c='darkgreen', s=12, label='predicted terminators', marker='x')

	plt.ylim(0,8.5)
	plt.xlim(0,12.5)
	plt.xlabel('Avg. RNA-Seq count')
	plt.ylabel('Max. Term-Seq count')
	plt.title('predicted classes')
	plt.legend()

	plt.savefig(outfileG1,dpi=300)

	# plt.show()














####################################################################################
####################################################################################

	plt.figure(figsize=(12,12),dpi=120)
	plt.suptitle('Own Decision boundary, PPV: ' + str(ppvG) + ', intercept (b): ' + str(b) + ', slope (m): ' + str(m))


	#plot training data (known classes), test data (known classes), test data all points (predicted classes)

	plt.subplot(221)

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='#E87CC6', alpha=0.8)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='#7A9441', alpha=0.5)


	plt.scatter(xvalsAllPoints[~maskII], yvalsAllPoints[~maskII], c='magenta', s=12, label='predicted negatives', marker='+')
	plt.scatter(xvalsAllPoints[maskII], yvalsAllPoints[maskII], c='#5AB91C', s=12, label='predicted terminators', marker='+')

	unique = list(set(classes))
	colors = ['darkmagenta', 'darkgreen']
	labelsPredicted = ['overlapping genes', 'overlapping terminators']

	for i, u in enumerate(unique):
	    xi = [Xtest[:,0][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]
	    yi = [Xtest[:,1][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, label=str(labelsPredicted[i]), marker='x')


	plt.ylim(0,8.5)
	plt.xlim(0,12.5)

	plt.xlabel('Avg. RNA-Seq count')
	plt.ylabel('Max. Term-Seq count')
	plt.title('predicted classes')

	plt.legend(bbox_to_anchor=(1.04,0.5), loc='center left')

################################################################################
	#plot training data (known classes), test data (known classes), test all points (predicted classes) overlapping test (red crosses)

	plt.subplot(223)

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='#E87CC6', alpha=0.8)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='#7A9441', alpha=0.5)



	unique = list(set(classes))
	colors = ['darkmagenta', 'darkgreen']
	labelsPredicted = ['overlapping terminators', 'overlapping genes']

	for i, u in enumerate(unique):
	    xi = [Xtest[:,0][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]
	    yi = [Xtest[:,1][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, marker='x')

	plt.scatter(npPredictedTerminatorsAllPointsOverlappingTest[:,0],npPredictedTerminatorsAllPointsOverlappingTest[:,1], c='red',s=11, label='overlapping', marker ='+')


	plt.ylim(0,8.5)
	plt.xlim(0,12.5)

	plt.xlabel('Avg. RNA-Seq count')
	plt.ylabel('Max. Term-Seq count')
	plt.title('predicted classes')

	plt.legend()

##############################################################################
	#plot training data (known classes), test data (known classes), test all points (predicted classes) NOT overlapping test (red crosses)

	plt.subplot(224)

	plt.plot(xII, yII, c='black',linewidth=1.0)
	plt.fill_between(xII, yII, color='#E87CC6', alpha=0.8)
	plt.fill_between(xII, yII, plt.ylim()[1],  color='#7A9441', alpha=0.5)


	unique = list(set(classes))
	colors = ['darkmagenta', 'darkgreen']
	labelsPredicted = ['overlapping terminators', 'overlapping genes']

	for i, u in enumerate(unique):
	    xi = [Xtest[:,0][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]
	    yi = [Xtest[:,1][j] for j  in range(len(Xtest[:,0])) if classes[j] == u]

	    plt.scatter(xi, yi, c=colors[i], s= 12, marker='x')

	plt.scatter(npPredictedTerminatorsAllPointsNOTOverlappingTest[:,0],npPredictedTerminatorsAllPointsNOTOverlappingTest[:,1], c='red',s=11, label='not overlapping', marker ='+')


	plt.ylim(0,8.5)
	plt.xlim(0,12.5)

	plt.xlabel('Avg. RNA-Seq count')
	plt.ylabel('Max. Term-Seq count')
	plt.title('predicted classes')

	plt.legend()

	plt.savefig(outfileG2,dpi=300)


