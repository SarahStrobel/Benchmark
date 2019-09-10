# scatter plots:

# 1)
# for each known terminator: take max 5'end position count of term-seq at these positions and max rna-seq coverage at point of highest term-seq
# or: take max end point count of term-seq at these positions and avg rna-seq coverage from this position to x nucs upstream

# 2)
# for whole genome: take max 5'end position count of term-seq and avg rna-seq count every x nucs; plot non-overlapping, overlapping genes and overlapping terminators


import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os.path
from operator import itemgetter


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
# methods for checking parsed file types

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def checkInt(v):
	v = int(v)
	if v < 0:
		raise argparse.ArgumentTypeError('positive Integer value expected')
	if isinstance(v, int):
		return v
	else:
		raise argparse.ArgumentTypeError('Integer value expected')

def checkIntList(v):
	v = v.split(' ')
	if len(v) != 2:
		raise argparse.ArgumentTypeError('two Integer values expected')

	v0 = int(v[0])
	v1 = int(v[1])

	if v0 < 0 or v1 < 0:
		raise argparse.ArgumentTypeError('positive Integer value expected')
	if v0 > v1:
		raise argparse.ArgumentTypeError('first number must be smaller than second number')
	if isinstance(v0, int) and isinstance(v1, int):
		return [v0,v1]
	else:
		raise argparse.ArgumentTypeError('Integer value expected')


def checkBedFormat(v):
	b = os.path.splitext(v)[1][1:].lower()
	if b != 'bed':
		raise argparse.ArgumentTypeError('bed format file type expected')
	else:
		return v

def checkBedgraphFormat(v):
	b = os.path.splitext(v)[1][1:].lower()
	if b != 'bedgraph':
		raise argparse.ArgumentTypeError('bedgraph format file type expected')
	else:
		return v

def checkBedgraphList(v):
	v = v.split(' ')
	bedgraphList = []
	for i in v:
		b = i
		b = os.path.splitext(b)[1][1:].lower()
		if b != 'bedgraph':
			raise argparse.ArgumentTypeError('bedgraph format file type expected')
		else:
			bedgraphList.append(i)
	return bedgraphList

def checkGffFormat(v):
	b = os.path.splitext(v)[1][1:].lower()
	if b != 'gff':
		raise argparse.ArgumentTypeError('gff format file type expected')
	else:
		return v


#######################################################################
#######################################################################
#checking if gene is big enough to chop off numberOfNucsToChopOffGenes, if not, substract 10 off numberOfNucsToChopOffGenes until it works

def subtractingNucs(geneLength, numberOfNucsToChopOffGenes):

	if geneLength <= 2*numberOfNucsToChopOffGenes:
		numberOfNucsToChopOffGenes = numberOfNucsToChopOffGenes - 10
		numberOfNucsToChopOffGenes = subtractingNucs(geneLength, numberOfNucsToChopOffGenes)
		return numberOfNucsToChopOffGenes

	else:
		return numberOfNucsToChopOffGenes



#######################################################################
#######################################################################
# parsing command line flags

parser = argparse.ArgumentParser(description= 'Get ratio of Maximum Term-Seq 5\' end nuc counts and avg./max RNA Seq coverage' + '\n'
								'Usage:' + '\t' + 'ratioTermSeq.py <options> -gff -bed -ts -rs')

#required files:
parser.add_argument('-gff', dest='gffFile', help='input of gene annotation file in gff format', type=checkGffFormat, required=True)
parser.add_argument('-bed', dest='terminatorBedFile', help='input of terminators in bed format', type=checkBedFormat, required=True)
parser.add_argument('-ts', dest='termSeqFileList', help='input of Term-Seq 5\' nucleotide frequency files in bedgraph format \
					(single file or space separated list of triplicates which will be averaged)', type=checkBedgraphList, required=True, nargs="+")
parser.add_argument('-rs', dest='rnaSeqFile', help='input of RNA-Seq coverage in bedgraph format',  type=checkBedgraphFormat, required=True)


#optional
parser.add_argument('-o', dest='outfile', help='output path and filename prefix, default: /ratioTermSeq', nargs='?', default="ratioTermSeq")
parser.add_argument('-max', dest='boolMax', help='calculate the maximum RNA-Seq coverage, default: False (calculate the mean RNA-Seq coverage)', type=str2bool, nargs='?', default=False)
parser.add_argument('-g', dest='numberOfNucsInGenome', help='use only the nucleotides from and to these positions in genome (space separated list, example: 1 500000),\
					 if the genome is smaller full genome will be used, default: 1 to 1,000,000', type=checkInt, default=[1,4215606],  nargs='+')
parser.add_argument('-c', dest='numberOfNucsToChopOffGenes', help='number of nucleotides to cut off both ends of genes, \
					if a gene is smaller size will be decremented by 10 until it works, default: 100', type=checkInt, nargs='?', default=100)
parser.add_argument('-s', dest='numberOfNucsToSplitInto', help='length of intervals to cut genome into, default: 50', type=checkInt, nargs='?', default=50)
parser.add_argument('-a', dest='numberOfNucsToAvg', help='length of region in RNA-Seq data to average/max over, default: 50', type=checkInt, nargs='?', default=50)



args = parser.parse_args()

gffFile = args.gffFile
terminatorBedFile = args.terminatorBedFile

termSeqFileList = args.termSeqFileList
if len(termSeqFileList) != 3:
	raise Exception('\n 3 replicates')

termSeqFile1 = ''
termSeqFile2 = ''
termSeqFile3 = ''

for i in range(len(termSeqFileList)):
	locals()["termSeqFile"+str(i+1)] = ''.join(termSeqFileList[i])


rnaSeqFile = args.rnaSeqFile
outfile = args.outfile
boolMax = args.boolMax

#######################################################################
#######################################################################

# numberOfNucsInGenome: only first x nucs
numberOfNucsInGenome = args.numberOfNucsInGenome
# numberOfNucsToChopOffGenes: chops off x nucs from both ends of genes
numberOfNucsToChopOffGenes = args.numberOfNucsToChopOffGenes
# numberOfNucsToSplitInto: divides the genome in parts of x nucs
numberOfNucsToSplitInto = args.numberOfNucsToSplitInto
#numberOfNucsToAvg: avg/max RNASeq count from the point of maxTScount to x nucs uptream
numberOfNucsToAvg = args.numberOfNucsToAvg
#average length of terminator
avgLengthOfTerminator = 0

outfile = outfile + '_' + str(numberOfNucsToAvg) + '_nucs'


print '\nspan of nucleotides in genome: \t\t\t\t' + str(numberOfNucsInGenome[0]) + '-' + str(numberOfNucsInGenome[1])
print 'number of nucleotides to cut off genes: \t\t' + str(numberOfNucsToChopOffGenes)
print 'length of intervals to cut the genome into: \t\t' + str(numberOfNucsToSplitInto)
print 'length of region in RNA-Seq to average/max over: \t' + str(numberOfNucsToAvg)
print 'name of outfile:\t\t\t\t\t' + str(outfile)
print 'calculate maximum RNA-Seq coverage:\t\t\t' + str(boolMax)

#######################################################################
#######################################################################
terminatorLengths = []
numberOfTerminators = 0
numberOfGenes = 0

terminatorCoords = []
geneCoords = []

maxTSandMeanRNATerminators = []
maxTSandMaxRNAWholeGenome = []
maxTSandMeanRNAWholeGenome = []
strandMaxTSandMeanRNAWholeGenome = []

maxTSstartAndCount = []
maxRSstartAndCount = []
maxTSstartAndCountWholeGenome = []
maxRSstartAndCountWholeGenome = []

TScountsOverlappingTerminators = []
TScountsOverlappingGenes = []

strandTScountsOverlappingTerminators = []
strandTScountsOverlappingGenes = []


#######################################################################
#######################################################################

outfileTSvsRS = outfile + '_TSvsRS'
outfileTScountsOverlappingTerminators = outfile + '_TScountsOverlappingTerminators'
otufileTScountsOverlappingGenes = outfile + '_TScountsOverlappingGenes'


#######################################################################
#######################################################################

with open(gffFile, 'r') as gff, open(terminatorBedFile, 'r') as bed, open(termSeqFile1, 'r') as ts, \
		open(termSeqFile2, 'r') as ts2, open(termSeqFile3, 'r') as ts3, open(rnaSeqFile, 'r') as rs:

	with open(outfileTSvsRS, 'w') as tsVsRs, open(outfileTScountsOverlappingTerminators, 'w') as tsOt, \
			open(otufileTScountsOverlappingGenes, 'w') as tsOg:

  
#######################################################################
#######################################################################
# # read only the span of "numberOfNucsInGenome = [from, to]" nucs of the Term-Seq (TS) and RNA-Seq (RS) files, if genome is smaller, use full genome

		# read all lines of the RNA-Seq file
		linesRS = rs.readlines()
		
		# read all lines of all 3 Term-Seq replicate files
		linesTS = ts.readlines()
		lines2TS = ts2.readlines()
		lines3TS = ts3.readlines()

		lengthGenome = len(linesRS)
		print 'lenght of genome:\t\t\t\t\t' + str(lengthGenome)

		# check if numbers in given span are in ascending order 		
		if numberOfNucsInGenome[0] > numberOfNucsInGenome[1]:
			raise Exception('\nfirst number for slicing genome greater than second number')
		# check if negative numbers
		if numberOfNucsInGenome[0] < 0 or numberOfNucsInGenome[1] < 0:
			raise Exception('\nnumbers for slicing genome must be greater than zero')

		# check if given span is smaller than length of the whole genome and slicing into given spans
		if numberOfNucsInGenome[0] < lengthGenome and numberOfNucsInGenome[1] <= lengthGenome:
			linesTS = linesTS[numberOfNucsInGenome[0]:numberOfNucsInGenome[1]]
			linesRS = linesRS[numberOfNucsInGenome[0]:numberOfNucsInGenome[1]]
			lines2TS = lines2TS[numberOfNucsInGenome[0]:numberOfNucsInGenome[1]]
			lines3TS = lines3TS[numberOfNucsInGenome[0]:numberOfNucsInGenome[1]]
		else:
			raise Exception('\nnucleotides to cut must be smaller than genome length')
			
		
		# lines of Term-Seq and RNA-Seq after slicing to given span
		linesTS2 = []
		lines2TS2 = []
		lines3TS2 = []

		linesRS2 = []


#######################################################################
#######################################################################
	# read gene annotation gff file and take start and end coordinates in the span of"numberOfNucsInGenome = [from, to]
	# and trim by "numberOfNucsToChopOffGenes" nucs on each side if the gene is big enough
	# if the gene is not big enough, decrease "numberOfNucsToChopOffGenes" recursively by 10 until it is

	# gff: 1-based, closed [start, end] 
	# from start to end, put all coords into a list to detect which points are overlapping known genes by binary search


		for line in gff:
			if '#' not in line:

				start = int(line.split()[3])
				end = int(line.split()[4])

				# check if start and end between the given span of nucleotides (numberOfNucsInGenome = [from, to])
				if start > numberOfNucsInGenome[0] and end < numberOfNucsInGenome[1]:

					# get start and end coords in given span, calculate length of gene/rna/cds
					start = int(line.split()[3])
					end = int(line.split()[4])

					geneLength = end - start

					if start > end:
						raise Exception('Gene annotation gff file: start smaller than end')

						
					# check if numberOfNucsToChopOff is small enough for each gene, if not, decrease recursively by 10 until it is
					newNumberOfNucsToChopOffGenes = subtractingNucs(geneLength, numberOfNucsToChopOffGenes)
					# adjust start and end coords
					start = start + newNumberOfNucsToChopOffGenes
					end = end - newNumberOfNucsToChopOffGenes



					if line.split()[2] == 'gene':
						# count number of genes in given span
						numberOfGenes = numberOfGenes + 1

						# put all gene coords between start and end into list:
						for i in range (start, end+1): #half open!
							geneCoords.append(i)


		
		#sort coord lists for binary search later
		geneCoords.sort()



#######################################################################
#######################################################################
	# read the terminator bed file line by line and take the start and end coordinates of each known terminator in the span of"numberOfNucsInGenome = [from, to]

	# bed: 0-based half open [start-1, end)
	# from startTerminator to endTerminator, put all coords into a list to detect which points are overlapping known terminators by binary search


		for line in bed:

			startTerminator = int(line.split()[1])+1 
			endTerminator = int(line.split()[2])+1

			# check if startTerminator and endTerminator between the given span of nucleotides (numberOfNucsInGenome = [from, to])
			if startTerminator > numberOfNucsInGenome[0] and endTerminator < numberOfNucsInGenome[1]:

				# count number of terminators, get start and end coords, calculate average length of terminators:
				numberOfTerminators = numberOfTerminators + 1

				startTerminator = int(line.split()[1]) + 1
				endTerminator = int(line.split()[2]) + 1
				lengthTerminator = endTerminator - startTerminator 

				terminatorLengths.append(lengthTerminator)

				avgLengthOfTerminator = np.mean(terminatorLengths)


				# put all coords between startTerminator and endTerminator into list:SS			
				for i in range (startTerminator, endTerminator):
					terminatorCoords.append(i)


				# for all Term-Seq replicates adjust linesTS2, lines2TS2 and lines3TS2 according to given span:
				linesTS2 = linesTS[(startTerminator -1 - numberOfNucsInGenome[0]) : (endTerminator -numberOfNucsInGenome[1])]
				lines2TS2 = lines2TS[(startTerminator -1 - numberOfNucsInGenome[0]) : (endTerminator -numberOfNucsInGenome[1])]
				lines3TS2 = lines3TS[(startTerminator -1 - numberOfNucsInGenome[0]) : (endTerminator -numberOfNucsInGenome[1])]


				###################################################
				#known terminators:
				#taking the maximum Term-seq count and its position between the coordinates of each known terminator 
		
				termSeqCountTerminators = {}
				termSeqStrandsTerminators = {}

				TSstart = 0
				# 3 replicates
				TScount = 0
				TScount2 = 0
				TScount3 = 0


				
				for i in range(len(linesTS2)):

					TSstart = int(linesTS2[i].split()[1])

					TScountStrand = []
					
					# check if TSstart between the given span of nucleotides (numberOfNucsInGenome = [from, to])
					if TSstart > numberOfNucsInGenome[0] and TSstart < numberOfNucsInGenome[1]:

						# get start position and Term-Seq count of each replicate and calculating the average/minimum
						TSstart = int(linesTS2[i].split()[1])

						TScount = int(linesTS2[i].split()[3])
						TScount2 = int(lines2TS2[i].split()[3])
						TScount3 = int(lines3TS2[i].split()[3])

						TSstrand = linesTS2[i].split()[-1]
						TSstrand2 = lines2TS2[i].split()[-1]
						TSstrand3 = lines3TS2[i].split()[-1]



						TScountStrand.append([TScount, TSstrand])
						TScountStrand.append([TScount2, TSstrand2])
						TScountStrand.append([TScount3, TSstrand3])

						maxTScount = max(TScountStrand,key=itemgetter(0))[0]
						maxTSstrand = max(TScountStrand,key=itemgetter(0))[1]


						# averaging/taking the minimum of the 3 replicates
						avgTScount = np.mean([TScount, TScount2, TScount3])
						minTScount = np.min([TScount, TScount2, TScount3])

						# print avgTScount

						# putting start position and average/minimum of Term-Seq counts in dictionary
						termSeqCountTerminators[TSstart] = avgTScount
						# termSeqCountTerminators[TSstart] = minTScount

						# 'avg' strand = strand of highest TScount
						termSeqStrandsTerminators[TSstart] = maxTSstrand

				# getting the maximum Term-Seq count and its position
				maxTScoord = max(termSeqCountTerminators.iterkeys(), key=lambda k: termSeqCountTerminators[k])
				maxTScount = termSeqCountTerminators[maxTScoord]
				maxTSstrand = termSeqStrandsTerminators[maxTScoord]


				# adjusting the position according to given span
				maxTScoord2 =  maxTScoord - numberOfNucsInGenome[0]


				# list of all coords and Term-Seq counts at a terminator's maxTScount (outliers)
				maxTSstartAndCount.append([maxTScoord,maxTScount])


				########################################################
				# known terminators:
				# taking the average/max RNA-seq coverage from the position of "maxTScoord" to "numberOfNucsToAvg" nucleotides downstream

				rnaSeqCountTerminators = []

				# slicing linesRS2 into spans of "numberOfNucsToAvg"
				linesRS2 = linesRS[maxTScoord2 - numberOfNucsToAvg : maxTScoord2]

				for i in range(len(linesRS2)):
					rnaSeqCountTerminators.append(int(linesRS2[i].split()[-1]))

				# either taking the mean or the max (depending on chosen flag) of the RNA-Seq counts in the span of "numberOfNucsToAvg"
				if boolMax == False:	
					RSCount = np.mean(rnaSeqCountTerminators)
				else:
					RSCount = np.max(rnaSeqCountTerminators)


				# maxTS and max RNA or avg RNA counts for terminators
				maxTSandMeanRNATerminators.append([maxTScount, RSCount])

				# list of all coords and RNA-Seq coverages at a terminator's maxTScount (outliers)
				maxRSstartAndCount.append([maxTScoord, RSCount, maxTScount])



		#sort terminator coord list for binary search later
		terminatorCoords.sort()


		# print out stats for terminators and genes
		print 'avg length of terminators:\t\t\t\t\t' + str(avgLengthOfTerminator)
		print 'number of terminators in the nucs ' + str(numberOfNucsInGenome) + ':\t\t' + str(numberOfTerminators) + ' (' + str(len(terminatorCoords)) + ' terminator nucleotides)'
		print 'number of genes in the nucs ' + str(numberOfNucsInGenome) + ':\t\t\t' + str(numberOfGenes) + ' (' + str(len(geneCoords)) + ' gene nucleotides)' 





#######################################################################
#######################################################################

	# maximum count for Term-Seq every "numberOfNucsToSplitInto" nucs and avgerage count for RNA-Seq from the position of maxTScoord to "numberOfNucsToAvg" nucs downstream

		# 3 replicates
		linesTS3 = []
		lines2TS3 = []
		lines3TS3 = []

		for i in range(0,(len(linesTS)), numberOfNucsToSplitInto):

			# slicing the Term-Seq count file in span of "numberOfNucsToSplitInto"
			# for example numberOfNucsToSplitInto = 50 --> first iteration: lines 0-50, second iteration lines 50-100, ...
			linesTS3 = linesTS[i: i+numberOfNucsToSplitInto]
			lines2TS3 = lines2TS[i: i+numberOfNucsToSplitInto]
			lines3TS3 = lines3TS[i: i+numberOfNucsToSplitInto]


			termSeqCountWholeGenome = {}
			TermSeqStrands = {}

			TSstart = 0
			TScount = 0
			TScount2 = 0
			TScount3 = 0


			for i in range(len(linesTS3)):

				TSstart = int(linesTS3[i].split()[1])

				TScountStrand = []

				# check if TSstart between the given span of nucleotides (numberOfNucsInGenome = [from, to])
				if TSstart > numberOfNucsInGenome[0] and TSstart < numberOfNucsInGenome[1]:

					# get start position and Term-Seq count of each replicate and calculating the average/minimum
					TSstart = int(linesTS3[i].split()[1])

					TScount = int(linesTS3[i].split()[3])
					TScount2 = int(lines2TS3[i].split()[3])
					TScount3 = int(lines3TS3[i].split()[3])

					TSstrand = linesTS3[i].split()[-1]
					TSstrand2 = lines2TS3[i].split()[-1]
					TSstrand3 = lines3TS3[i].split()[-1]

					TScountStrand.append([TScount, TSstrand])
					TScountStrand.append([TScount2, TSstrand2])
					TScountStrand.append([TScount3, TSstrand3])

					maxTScount = max(TScountStrand,key=itemgetter(0))[0]
					maxTSstrand = max(TScountStrand,key=itemgetter(0))[1]


					# averaging/taking the minimum of the 3 replicates
					avgTScount = np.mean([TScount, TScount2, TScount3])
					minTScount = np.min([TScount, TScount2, TScount3])

					# putting start position and average/minimum of Term-Seq counts in dictionary
					termSeqCountWholeGenome[TSstart] = avgTScount
					TermSeqStrands[TSstart] = maxTSstrand
			

			maxTScoord = max(termSeqCountWholeGenome.iterkeys(), key=lambda k: termSeqCountWholeGenome[k])
			maxTScount = termSeqCountWholeGenome[maxTScoord]
			maxTSstrand = TermSeqStrands[maxTScoord]


			# adjusting the position according to given span
			maxTScoord2 =  maxTScoord - numberOfNucsInGenome[0]

			# for getting all coords and counts at a terminator's maxTScount
			maxTSstartAndCountWholeGenome.append([maxTScoord,maxTScount])



		 	######################################################
			# taking the RNA-seq max/average coverage from the position of maxTScoord to "numberOfNucsToAvg" nucs downstream
			
			rnaSeqCountWholeGenome = []

			# if the maxTScoord is too small, subtracting could lead to negative values
			if maxTScoord2 - numberOfNucsToAvg < 0:
				linesRS3 = linesRS[0 : maxTScoord2]

			else:
				linesRS3 = linesRS[maxTScoord2 - numberOfNucsToAvg : maxTScoord2]


			for i in range(len(linesRS3)):
				rnaSeqCountWholeGenome.append(int(linesRS3[i].split()[-1]))


			# either taking the mean or the max (depending on chosen flag) of the RNA-Seq counts in the span of "numberOfNucsToAvg"
			if boolMax == False:
				RSCountWholeGenome = np.mean(rnaSeqCountWholeGenome)
			else:
				RSCountWholeGenome = np.max(rnaSeqCountWholeGenome)

			

			# maxTS and max RNA or avg RNA counts for whole genome
			maxTSandMeanRNAWholeGenome.append([maxTScoord, maxTScount, RSCountWholeGenome])
			strandMaxTSandMeanRNAWholeGenome.append([maxTScoord, maxTScount, RSCountWholeGenome,maxTSstrand])

			# for getting all maxTScoords and avgRScounts at a terminator's maxTScount (outliers)
			maxRSstartAndCountWholeGenome.append([maxTScoord, RSCountWholeGenome, maxTScount])



	######################################################
	# binary search for checking if points overlap known terminators / predicted terminators / known genes/CDS/RNA


			resultTerminators = binarySearch(terminatorCoords, 0, len(terminatorCoords)-1, maxTScoord)

			if resultTerminators != -1:
				TScountsOverlappingTerminators.append([maxTScoord, maxTScount,RSCountWholeGenome])
				strandTScountsOverlappingTerminators.append([maxTScoord, maxTScount,RSCountWholeGenome, maxTSstrand])

			
			resultsGenes = binarySearch(geneCoords, 0, len(geneCoords)-1, maxTScoord)

			if resultsGenes != -1:
				TScountsOverlappingGenes.append([maxTScoord, maxTScount, RSCountWholeGenome])
				strandTScountsOverlappingGenes.append([maxTScoord, maxTScount, RSCountWholeGenome, maxTSstrand])


		# if points overlap predicted terminators or points overlap genes:				
		TScountsOverlappingTerminatorsPLUSTScountsOverlappingGenes = TScountsOverlappingTerminators + TScountsOverlappingGenes
		strandTScountsOverlappingTerminatorsPLUSTScountsOverlappingGenes = strandTScountsOverlappingTerminators + strandTScountsOverlappingGenes



#######################################################################
#######################################################################
		# print stats for overlaps

		print '\nTS overl. known Terminators: ' + str(len(TScountsOverlappingTerminators))
		print 'TS overl. Genes: ' + str(len(TScountsOverlappingGenes))
		print '\n'

#######################################################################
#######################################################################
	# write to files with coords and strand 


		for list1 in strandMaxTSandMeanRNAWholeGenome:
			tsVsRs.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n') # with coords and strand


		for list1 in strandTScountsOverlappingTerminators:
			tsOt.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n')


		for list1 in strandTScountsOverlappingGenes:
			tsOg.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n')



#######################################################################
#######################################################################
# plotting: 

	npMaxTSandMeanRNATerminators = np.array(maxTSandMeanRNATerminators)

	plt.figure(1, figsize =(12,6), dpi=100)

	# at known terminators
	sub1 = plt.subplot(121) 

	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(1, 100000)
	plt.ylim(1, 100000)
	plt.grid(True)
	# adding +1 to all values --> log not - infinity
	plt.scatter(npMaxTSandMeanRNATerminators[:,1]+1, npMaxTSandMeanRNATerminators[:,0]+1, s = 15, color = 'red')
	plt.title("Term-Seq vs. RNA-Seq at known terminators \n(avg of " + str(numberOfNucsToAvg) + " nucleotide intervals)")
	plt.ylabel("Max. Term-Seq count")
	plt.xlabel("Avg. RNA-Seq count")

	sub1.set_axisbelow(True)


	npMaxTSandMeanRNAWholeGenome = np.array(maxTSandMeanRNAWholeGenome)
	npTSCountsOverlappingTerminators = np.array(TScountsOverlappingTerminators)
	npTScountsOverlappingGenes = np.array(TScountsOverlappingGenes)


	sub2 = plt.subplot(122) 

	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(1, 100000)
	plt.ylim(1, 100000)

	plt.scatter(npMaxTSandMeanRNAWholeGenome[:,2]+1, npMaxTSandMeanRNAWholeGenome[:,1]+1, s = 15, label = "non-overlapping")
	plt.scatter(npTScountsOverlappingGenes[:,2]+1, npTScountsOverlappingGenes[:,1]+1, s = 15, color = 'plum', label = "overlapping genes")
	plt.scatter(npTSCountsOverlappingTerminators[:,2]+1, npTSCountsOverlappingTerminators[:,1]+1, s = 15, color = 'red', label = "overlapping known terminators")

	plt.legend()
	plt.grid(True)
	plt.title("Term-Seq vs. RNA-Seq \n(avg of " + str(numberOfNucsToAvg) + " nucleotide intervals)")
	plt.ylabel("Max. Term-Seq count")
	plt.xlabel("Avg. RNA-Seq count")

	sub2.set_axisbelow(True)


	plt.savefig(outfile, dpi=300)


	plt.close()

#######################################################################
#######################################################################



gff.close()
bed.close()
ts.close()
rs.close()


tsVsRs.close()

tsOt.close()
tsOg.close()
