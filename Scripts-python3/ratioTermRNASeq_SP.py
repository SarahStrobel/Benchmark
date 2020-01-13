# makes scatter plots of avg term-seq 5' end nuc counts vs. max rna-seq counts
# writes files for all points, points overlapping known terminators and points overlapping genes (B.subtilis)
# writes files for all points and points overlapping genes (E.faecalis + L.monocytogenes)
# writes files for RNIE making 'artificial terminators' of length 120 (E.faecalis + L.monocytogenes)
# prints statistics


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import sys
import os.path
from operator import itemgetter
import glob
import math


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
	if not (b == 'gff' or b == 'gff3'):
		raise argparse.ArgumentTypeError('gff or gff3 format file type expected')
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
								'Usage:' + '\t' + 'ratioTermSeq.py <options> -gff -bed -ts -rs' +'\n'
								'optional:' + '\t' + '-o -max -g -c -s - a -l')

#required files:
parser.add_argument('-gff', dest='gffFile', help='input of gene annotation file in gff format', type=checkGffFormat, required=True)
parser.add_argument('-bed', dest='terminatorBedFile', help='input of terminators in bed format', type=checkBedFormat, required=True)
parser.add_argument('-ts', dest='termSeqFileList', help='input of Term-Seq 5\' nucleotide frequency files in bedgraph format \
					(single file or space separated list of triplicates which will be averaged)', type=checkBedgraphList, required=True, nargs="+")
parser.add_argument('-rs', dest='rnaSeqFile', help='input of RNA-Seq coverage in bedgraph format',  type=checkBedgraphFormat, required=True)


#optional
parser.add_argument('-o', dest='outfile', help='output path and filename prefix, default: /ratioTermSeq', nargs='?', default="ratioTermSeq")
parser.add_argument('-max', dest='boolMax', help='calculate the maximum RNA-Seq coverage, default: False (calculate the mean RNA-Seq coverage)', type=str2bool, nargs='?', default=False)
parser.add_argument('-g', dest='numberOfNucsInGenome', help='use only the nucleotides from and to these positions in genome (zero based, space separated list, example: 0 500000),\
					 if the genome is smaller full genome will be used, default: 1 to 1,000,000', type=checkInt, default=[],  nargs='+')
parser.add_argument('-c', dest='numberOfNucsToChopOffGenes', help='number of nucleotides to cut off both ends of genes, \
					if a gene is smaller size will be decremented by 10 until it works, default: 0', type=checkInt, nargs='?', default=0)
parser.add_argument('-s', dest='numberOfNucsToSplitInto', help='length of intervals to cut genome into, default: 50', type=checkInt, nargs='?', default=50)
parser.add_argument('-a', dest='numberOfNucsToAvg', help='length of region in RNA-Seq data to average/max over, default: 50', type=checkInt, nargs='?', default=50)
parser.add_argument('-l', dest='lengthTerminator', help='length of artificial terminator, default:100', type=checkInt, nargs='?', default=100)


args = parser.parse_args()

gffFile = args.gffFile
terminatorBedFile = args.terminatorBedFile

termSeqFileList = args.termSeqFileList
if not (len(termSeqFileList) == 3 or len(termSeqFileList) == 4):
	raise Exception('\n 3 or 4 replicates')


termSeqFileList = [item for sublist in termSeqFileList for item in sublist] #flatten list
termSeqPath = os.path.dirname(termSeqFileList[1]) + '/'


rnaSeqFile = args.rnaSeqFile
rnaSeqPath = os.path.dirname(rnaSeqFile) + '/'


outfile = args.outfile
boolMax = args.boolMax

lengthTerminator = args.lengthTerminator

l1 = lengthTerminator * 0.16666666666666664
l2 = lengthTerminator - l1

l1 =  int(math.ceil(l1))
l2 = int(math.floor(l2))

# print l1
# print l2

# ####################################################################################################################
bsubTermSeqFiles = []

experiments = (termSeqPath + '*ERX1304415*.bedgraph', termSeqPath + '*ERX1320300*.bedgraph', termSeqPath + '*ERX1320301*.bedgraph')
for experiment in experiments:
	bsubTermSeqFiles.append(glob.glob(experiment))

bsubTermSeqFiles = [item for sublist in bsubTermSeqFiles for item in sublist]

bsubRNASeqFile = rnaSeqPath + "sorted_filtered_trimmed_ERX1320302_Bacillus_subtilis_genomeCoverageBed.bedgraph"


lmonTermSeqFiles = []

experiments = (termSeqPath + '*ERR1248436*.bedgraph', termSeqPath + '*ERR1248437*.bedgraph', termSeqPath + '*ERR1248438*.bedgraph')
for experiment in experiments:
	lmonTermSeqFiles.append(glob.glob(experiment))

lmonTermSeqFiles = [item for sublist in lmonTermSeqFiles for item in sublist]



efaecTermSeqFiles = []

for bedgraphFile in glob.glob(termSeqPath + '*_chromosome.bedgraph'):
	efaecTermSeqFiles.append(str(bedgraphFile))

for bedgraphFile in glob.glob(termSeqPath + '*_plasmid1.bedgraph'):
	efaecTermSeqFiles.append(str(bedgraphFile))

for bedgraphFile in glob.glob(termSeqPath + '*_plasmid2.bedgraph'):
	efaecTermSeqFiles.append(str(bedgraphFile))

for bedgraphFile in glob.glob(termSeqPath + '*_plasmid3.bedgraph'):
	efaecTermSeqFiles.append(str(bedgraphFile))



efaecRNASeqFiles = []
for file in glob.glob(rnaSeqPath + 'sorted_filtered_trimmed_ERR1248404*.bedgraph'):
	efaecRNASeqFiles.append(str(file))




spneuTermSeqFiles = []

experiments = (termSeqPath + '*SRR7160964*.bedgraph', termSeqPath + '*SRR7160965*.bedgraph', termSeqPath + '*SRR7160966*.bedgraph', termSeqPath + '*SRR7160967*.bedgraph')
for experiment in experiments:
	spneuTermSeqFiles.append(glob.glob(experiment))

spneuTermSeqFiles = [item for sublist in spneuTermSeqFiles for item in sublist]

#########################################################################################################################

organism = ''
chrom = ''
plasmid = ''

if 'Bacillus' in rnaSeqFile:
	organism = 'B.subtilis'
	chrom = 'NC_000964.3'
	plasmid = 'chromosome'

if 'Listeria' in rnaSeqFile:
	organism = 'L.monocytogenes'
	chrom = "NC_003210.1"
	plasmid = 'chromosome'


if 'Enterococcus' in rnaSeqFile:
	organism = 'E.faecalis'
	if '_chromosome' in rnaSeqFile:
		chrom = "NC_004668.1"
		plasmid = 'chromosome'
		print(plasmid)
	if 'plasmid1' in rnaSeqFile:
		chrom = "NC_004669.1"
		plasmid = 'plasmid1'
		print(plasmid)
	if 'plasmid2' in rnaSeqFile:
		chrom = "NC_004671.1"
		plasmid = 'plasmid2'
		print(plasmid)
	if 'plasmid3' in rnaSeqFile:
		chrom = "NC_004670.1"
		plasmid = 'plasmid3'
		print(plasmid)

if 'Streptococcus' in rnaSeqFile:
	organism = 'S.pneumoniae'
	chrom = "NC_003028.3"
	plasmid = 'chromosome'


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
#average/max length of terminators
avgLengthOfTerminator = 0
maxLengthOfTerminator = 0

outfile = outfile + '_' + str(numberOfNucsToAvg) + '_nucs'


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
outfileTScountsOverlappingGenes = outfile + '_TScountsOverlappingGenes'

RNIEoutfileTSvsRS = outfile + 'RNIE_TSvsRS.bed'
RNIEoutfileTScountsOverlappingTerminators = outfile + 'RNIE_TScountsOverlappingTerminators.bed'
RNIEoutfileTScountsOverlappingGenes = outfile + 'RNIE_TScountsOverlappingGenes.bed'

#######################################################################
#######################################################################
#count unique reads in b.sub for normalization against unique reads in other organisms

uniqueReadsRSBS = 0
uniqueReadsTSCombinedBS = 0

uniqueReadsRS = 0
uniqueReadsTSCombined = 0


bsubFilesTS = [open(i, 'r') for i in bsubTermSeqFiles]
bsubFileRS = open(bsubRNASeqFile, 'r')

for file in bsubFilesTS:
	headerTS = file.readline()

	uniqueReads = int(headerTS.split('"')[1])
	uniqueReadsTSCombinedBS += uniqueReads

headerRS = bsubFileRS.readline()
uniqueReadsRSBS = int(headerRS.split('"')[1])



if organism == 'B.subtilis':
	uniqueReadsTSCombined = uniqueReadsTSCombinedBS
	uniqueReadsRS = uniqueReadsRSBS


if organism == 'L.monocytogenes':
	lmonFilesTS = [open(i, 'r') for i in lmonTermSeqFiles]
	lmonFileRS = open(rnaSeqFile, 'r')

	for file in lmonFilesTS:
		headerTS = file.readline()

		uniqueReads = int(headerTS.split('"')[1])
		uniqueReadsTSCombined += uniqueReads

	headerRS = lmonFileRS.readline()
	uniqueReadsRS = int(headerRS.split('"')[1])


if organism == 'E.faecalis':
	efaecFilesTS = [open(i, 'r') for i in efaecTermSeqFiles]
	efaecFilesRS = [open(i, 'r') for i in efaecRNASeqFiles]

	for file in efaecFilesTS:
		headerTS = file.readline()

		uniqueReads = int(headerTS.split('"')[1])
		uniqueReadsTSCombined += uniqueReads

	for file in efaecFilesRS:
		headerRS = file.readline()

		uniqueReadsRS += int(headerRS.split('"')[1])

if organism == 'S.pneumoniae':
	print(organism)
	spneuFilesTS = [open(i, 'r') for i in spneuTermSeqFiles]
	spneuFileRS = open(rnaSeqFile, 'r')

	for file in spneuFilesTS:
		headerTS = file.readline()

		uniqueReads = int(headerTS.split('"')[1])
		uniqueReadsTSCombined += uniqueReads

	headerRS = spneuFileRS.readline()
	uniqueReadsRS = int(headerRS.split('"')[1])
	
########################################################################
# open Term-Seq replicates, gene annotation file, terminator file and RNA-Seq file	
# read only the span of "numberOfNucsInGenome = [from, to]" nucs of the Term-Seq (TS) and RNA-Seq (RS) files

termSeqFiles = [open(i, 'r') for i in termSeqFileList]

with open(gffFile, 'r') as gff, open(terminatorBedFile, 'r') as bed, \
	open(rnaSeqFile, 'r') as rs: \
	
	#skipping header
	rs.readline()

	linesRS = rs.readlines()
	lengthGenome = len(linesRS)

	# if there is no span for -g option --> default to whole genome
	if not numberOfNucsInGenome:
		numberOfNucsInGenome.append(0)
		numberOfNucsInGenome.append(lengthGenome)

	# check if numbers in given span are in ascending order 		
	if numberOfNucsInGenome[0] > numberOfNucsInGenome[1]:
		raise Exception('\nfirst number for slicing genome greater than second number')
	# check if negative numbers
	if numberOfNucsInGenome[0] < 0 or numberOfNucsInGenome[1] < 0:
		raise Exception('\nnumbers for slicing genome must be greater than zero')
	# check if given span is smaller than length of the whole genome and slicing into given spans
	if numberOfNucsInGenome[0] < lengthGenome and numberOfNucsInGenome[1] <= lengthGenome:
		linesRS = linesRS[numberOfNucsInGenome[0]:numberOfNucsInGenome[1]]
	else:
		raise Exception('\nnucleotides to cut must be smaller than genome length')		

	print('\norganism:\t\t\t\t\t\t' + organism)
	print('unique reads in all Term-Seq files:\t\t\t' + str(uniqueReadsTSCombined))
	print('unique reads in RNA-Seq files\t\t\t\t' + str(uniqueReadsRS))
	print('span of nucleotides in genome: \t\t\t\t' + str(numberOfNucsInGenome[0]) + '-' + str(numberOfNucsInGenome[1]))
	print('number of nucleotides to cut off genes: \t\t' + str(numberOfNucsToChopOffGenes))
	print('length of intervals to cut the genome into: \t\t' + str(numberOfNucsToSplitInto))
	print('length of region in RNA-Seq to average/max over: \t' + str(numberOfNucsToAvg))
	print('path of outfile:\t\t\t\t\t' + str(outfile))
	print('calculate maximum RNA-Seq coverage:\t\t\t' + str(boolMax))
	print('lenght of genome:\t\t\t\t\t' + str(lengthGenome))


	linesTSall = []

	for ts in termSeqFiles:
		# skip header
		ts.readline()

		linesTS = ts.readlines()
		
		# check if given span is smaller than length of the whole genome and slicing into given span
		if numberOfNucsInGenome[0] < lengthGenome and numberOfNucsInGenome[1] <= lengthGenome:
			linesTS = linesTS[numberOfNucsInGenome[0]:numberOfNucsInGenome[1]]
			linesTSall.append(linesTS)
		else:
			raise Exception('\nnucleotides to cut must be smaller than genome length')
			
	linesTS = linesTSall[0]
	lines2TS = linesTSall[1]
	lines3TS = linesTSall[2]
	lines4TS = linesTSall[3]


#######################################################################
#######################################################################
# read gene annotation gff file and take start and end coordinates in the span of"numberOfNucsInGenome = [from, to]
# trim by "numberOfNucsToChopOffGenes" nucs on each side if the gene is big enough
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

	# lines of Term-Seq and RNA-Seq after slicing into given spans (default 50)
	linesRS2 = []
	linesTS2 = []
	lines2TS2 = []
	lines3TS2 = []
	lines4TS2 = []

	if organism == 'B.subtilis':
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
				maxLengthOfTerminator = np.max(terminatorLengths)


				# put all coords between startTerminator and endTerminator into list:SS			
				for i in range (startTerminator, endTerminator):
					terminatorCoords.append(i)


				# for all Term-Seq replicates adjust linesTS2, lines2TS2 and lines3TS2 according to given span:
				linesTS2 = linesTS[(startTerminator -1 - numberOfNucsInGenome[0]) : (endTerminator -numberOfNucsInGenome[1])]
				lines2TS2 = lines2TS[(startTerminator -1 - numberOfNucsInGenome[0]) : (endTerminator -numberOfNucsInGenome[1])]
				lines3TS2 = lines3TS[(startTerminator -1 - numberOfNucsInGenome[0]) : (endTerminator -numberOfNucsInGenome[1])]
				lines4TS2 = lines4TS[(startTerminator -1 - numberOfNucsInGenome[0]) : (endTerminator -numberOfNucsInGenome[1])]

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
				TScount4 = 0

				
				for i in range(len(linesTS2)):

					TSstart = int(linesTS2[i].split()[1])

					TScountStrand = []
					
					# check if TSstart between the given span of nucleotides (numberOfNucsInGenome = [from, to])
					if TSstart > numberOfNucsInGenome[0] and TSstart < numberOfNucsInGenome[1]:

						# get start position and Term-Seq count of each replicate and calculating the average/minimum
						TSstart = int(linesTS2[i].split()[1])

						#normalized: counts*number of unique reads in B.sub / number of unique reads in EF(or LM)

						TScount = int(linesTS2[i].split()[3]) * uniqueReadsTSCombinedBS / uniqueReadsTSCombined
						TScount2 = int(lines2TS2[i].split()[3]) * uniqueReadsTSCombinedBS / uniqueReadsTSCombined
						TScount3 = int(lines3TS2[i].split()[3])	* uniqueReadsTSCombinedBS / uniqueReadsTSCombined
						TScount4 = int(lines4TS2[i].split()[3])	* uniqueReadsTSCombinedBS / uniqueReadsTSCombined

						TSstrand = linesTS2[i].split()[-1]
						TSstrand2 = lines2TS2[i].split()[-1]
						TSstrand3 = lines3TS2[i].split()[-1]
						TSstrand4 = lines4TS2[i].split()[-1]


						TScountStrand.append([TScount, TSstrand])
						TScountStrand.append([TScount2, TSstrand2])
						TScountStrand.append([TScount3, TSstrand3])
						TScountStrand.append([TScount4, TSstrand4])

						maxTScount = max(TScountStrand,key=itemgetter(0))[0]
						maxTSstrand = max(TScountStrand,key=itemgetter(0))[1]


						# averaging/taking the minimum of the 3 replicates
						# avgTScount = np.mean([TScount, TScount2, TScount3])
						# minTScount = np.min([TScount, TScount2, TScount3])
						avgTScount = np.mean([TScount, TScount2, TScount3, TScount4])

						# putting start position and average/minimum of Term-Seq counts in dictionary
						termSeqCountTerminators[TSstart] = avgTScount
						# termSeqCountTerminators[TSstart] = minTScount

						# 'avg' strand = strand of highest TScount
						termSeqStrandsTerminators[TSstart] = maxTSstrand

				# getting the maximum Term-Seq count and its position
				maxTScoord = max(iter(termSeqCountTerminators.keys()), key=lambda k: termSeqCountTerminators[k])
				maxTScount = termSeqCountTerminators[maxTScoord]
				maxTSstrand = termSeqStrandsTerminators[maxTScoord]


				# adjusting the position according to given span
				maxTScoord2 =  maxTScoord - numberOfNucsInGenome[0]


				########################################################
				# known terminators:
				# taking the average/max RNA-seq coverage from the position of "maxTScoord" to "numberOfNucsToAvg" nucleotides downstream

				rnaSeqCountTerminators = []
				rnaSeqCountTerminators1 = []
				rnaSeqCountTerminators2 = []

				RSCount = 0.0

				#slicing linesRS into spans of "numberOfNucsToAvg" according to strand
				if maxTSstrand == '-':
					linesRS2 = linesRS[maxTScoord2 - numberOfNucsToAvg : maxTScoord2]
				if maxTSstrand == '+':
					linesRS2 = linesRS[maxTScoord2 : maxTScoord2 + numberOfNucsToAvg]

				# if there is no strand info, go 50 downstream/upstream and choose minimum of the two 
				else:
					linesRS21 = linesRS[maxTScoord2 - numberOfNucsToAvg : maxTScoord2]
					linesRS22 = linesRS[maxTScoord2 : maxTScoord2 + numberOfNucsToAvg]

					for i in range(len(linesRS21)):
						rnaSeqCountTerminators1.append(int(linesRS21[i].split()[-1]))
					for i in range(len(linesRS22)):
						rnaSeqCountTerminators2.append(int(linesRS22[i].split()[-1]))

					if boolMax == False:	
						RSCount = np.min([np.mean(rnaSeqCountTerminators1), np.mean(rnaSeqCountTerminators2)]) * uniqueReadsRSBS / uniqueReadsRS
					else:
						RSCount = np.min([np.max(rnaSeqCountTerminators1), np.max(rnaSeqCountTerminators2)]) * uniqueReadsRSBS / uniqueReadsRS




				# append RNASeq counts of all positions between maxTScoord and "numberOfNucsToAvg" downstream/upstream
				for i in range(len(linesRS2)):
					rnaSeqCountTerminators.append(int(linesRS2[i].split()[-1]))


				if boolMax == False:	
					RSCount = np.mean(rnaSeqCountTerminators) * uniqueReadsRSBS / uniqueReadsRS
				else:
					RSCount = np.max(rnaSeqCountTerminators) * uniqueReadsRSBS / uniqueReadsRS


				# maxTS and max RNA or avg RNA counts for terminators
				maxTSandMeanRNATerminators.append([maxTScount, RSCount])



		#sort terminator coord list for binary search later
		terminatorCoords.sort()


		# print out stats for terminators and genes
		print('max lenght of terminators:\t\t\t\t' + str(maxLengthOfTerminator))
		print('avg length of terminators:\t\t\t\t' + str(avgLengthOfTerminator))
		print('number of terminators in the nucs ' + str(numberOfNucsInGenome) + ':\t\t' + str(numberOfTerminators) + ' (' + str(len(terminatorCoords)) + ' terminator nucleotides)')


	print('number of genes in the nucs ' + str(numberOfNucsInGenome) + ':\t\t' + str(numberOfGenes) + ' (' + str(len(geneCoords)) + ' gene nucleotides)') 





#######################################################################
#######################################################################

# maximum count for Term-Seq every "numberOfNucsToSplitInto" nucs and avgerage count for RNA-Seq from the position of maxTScoord to "numberOfNucsToAvg" nucs downstream

	# 3 replicates
	linesTS3 = []
	lines2TS3 = []
	lines3TS3 = []
	linees4TS3 = []


	linesRS3 = []
	linesRS31 = []
	linesRS32 = []
	linesR432 = []

	for i in range(0,(len(linesTS)), numberOfNucsToSplitInto):

		# slicing the Term-Seq count file in span of "numberOfNucsToSplitInto"
		# for example numberOfNucsToSplitInto = 50
		linesTS3 = linesTS[i: i+numberOfNucsToSplitInto]
		lines2TS3 = lines2TS[i: i+numberOfNucsToSplitInto]
		lines3TS3 = lines3TS[i: i+numberOfNucsToSplitInto]
		lines4TS3 = lines4TS[i: i+numberOfNucsToSplitInto]


		termSeqCountWholeGenome = {}
		TermSeqStrands = {}

		TSstart = 0
		TScount = 0
		TScount2 = 0
		TScount3 = 0
		TScount4 = 0


		for i in range(len(linesTS3)):

			TSstart = int(linesTS3[i].split()[1])

			TScountStrand = []

			# check if TSstart between the given span of nucleotides (numberOfNucsInGenome = [from, to])
			if TSstart > numberOfNucsInGenome[0] and TSstart < numberOfNucsInGenome[1]:

				# get start position and Term-Seq count of each replicate and calculating the average/minimum
				TSstart = int(linesTS3[i].split()[1])

				TScount = int(linesTS3[i].split()[3]) * uniqueReadsTSCombinedBS / uniqueReadsTSCombined
				TScount2 = int(lines2TS3[i].split()[3]) * uniqueReadsTSCombinedBS / uniqueReadsTSCombined
				TScount3 = int(lines3TS3[i].split()[3]) * uniqueReadsTSCombinedBS / uniqueReadsTSCombined
				TScount4 = int(lines4TS3[i].split()[3]) * uniqueReadsTSCombinedBS / uniqueReadsTSCombined

				TSstrand = linesTS3[i].split()[-1]
				TSstrand2 = lines2TS3[i].split()[-1]
				TSstrand3 = lines3TS3[i].split()[-1]
				TSstrand4 = lines4TS3[i].split()[-1]

				TScountStrand.append([TScount, TSstrand])
				TScountStrand.append([TScount2, TSstrand2])
				TScountStrand.append([TScount3, TSstrand3])
				TScountStrand.append([TScount4, TSstrand4])

				maxTScount = max(TScountStrand,key=itemgetter(0))[0]
				maxTSstrand = max(TScountStrand,key=itemgetter(0))[1]


				# averaging/taking the minimum of the 3 replicates
				# avgTScount = np.mean([TScount, TScount2, TScount3])
				# minTScount = np.min([TScount, TScount2, TScount3])
				avgTScount = np.mean([TScount, TScount2, TScount3, TScount4])

				# putting start position and average/minimum of Term-Seq counts in dictionary
				termSeqCountWholeGenome[TSstart] = avgTScount
				TermSeqStrands[TSstart] = maxTSstrand
		

		maxTScoord = max(iter(termSeqCountWholeGenome.keys()), key=lambda k: termSeqCountWholeGenome[k])
		maxTScount = termSeqCountWholeGenome[maxTScoord]
		maxTSstrand = TermSeqStrands[maxTScoord]


		# adjusting the position according to given span
		maxTScoord2 =  maxTScoord - numberOfNucsInGenome[0]

		# for getting all coords and counts at a terminator's maxTScount
		maxTSstartAndCountWholeGenome.append([maxTScoord,maxTScount])



	 	######################################################
		# taking the RNA-seq max/average coverage from the position of maxTScoord to "numberOfNucsToAvg" nucs downstream/upstream
		
		rnaSeqCountWholeGenome = []

		rnaSeqCountWholeGenome1 = []
		rnaSeqCountWholeGenome2 = []

		RSCountWholeGenome = 0.0


		# slicing linesRS into spans of "numberOfNucsToAvg" according to strand
		if maxTSstrand == '-':
			# if the maxTScoord is too small, subtracting could lead to negative values
			if maxTScoord2 - numberOfNucsToAvg < 0:
				linesRS3 = linesRS[0 : maxTScoord2]
			else:
				linesRS3 = linesRS[maxTScoord2 - numberOfNucsToAvg : maxTScoord2]

			for i in range(len(linesRS3)):
				rnaSeqCountWholeGenome.append(int(linesRS3[i].split()[-1]))
			# either taking the mean or the max (depending on chosen flag) of the RNA-Seq counts in the span of "numberOfNucsToAvg"	
			if boolMax == False:
				RSCountWholeGenome = np.mean(rnaSeqCountWholeGenome) * uniqueReadsRSBS / uniqueReadsRS
			else:
				RSCountWholeGenome = np.max(rnaSeqCountWholeGenome) * uniqueReadsRSBS / uniqueReadsRS

		elif maxTSstrand == '+':
			# if the maxTScoord is too big, adding could lead to values bigger than genome length
			if maxTScoord2 + numberOfNucsToAvg > lengthGenome:
				linesRS3 = linesRS[maxTScoord2 : lengthGenome]
			else:	
				linesRS3 = linesRS[maxTScoord2 : maxTScoord2 + numberOfNucsToAvg]

			for i in range(len(linesRS3)):
				rnaSeqCountWholeGenome.append(int(linesRS3[i].split()[-1]))
			# either taking the mean or the max (depending on chosen flag) of the RNA-Seq counts in the span of "numberOfNucsToAvg"	
			if boolMax == False:
				RSCountWholeGenome = np.mean(rnaSeqCountWholeGenome) * uniqueReadsRSBS / uniqueReadsRS
			else:
				RSCountWholeGenome = np.max(rnaSeqCountWholeGenome) * uniqueReadsRSBS / uniqueReadsRS

		# if there is no strand info, go 50 downstream/upstream and choose minimum of the two  
		else:
			if maxTScoord2 - numberOfNucsToAvg < 0:
				linesRS31 = linesRS[0 : maxTScoord2]				
			else:
				linesRS31 = linesRS[maxTScoord2 - numberOfNucsToAvg : maxTScoord2]
			if maxTScoord2 + numberOfNucsToAvg > lengthGenome:
				linesRS32 = linesRS[maxTScoord2 : lengthGenome]
			else:
				linesRS32 = linesRS[maxTScoord2 : numberOfNucsToAvg + maxTScoord2]

			
			for i in range(len(linesRS31)):
				rnaSeqCountWholeGenome1.append(int(linesRS31[i].split()[-1])) 
			for i in range(len(linesRS32)):
				rnaSeqCountWholeGenome2.append(int(linesRS32[i].split()[-1])) 

			if boolMax == False:	
				RSCountWholeGenome = np.min([np.mean(rnaSeqCountWholeGenome1), np.mean(rnaSeqCountWholeGenome2)])* uniqueReadsRSBS / uniqueReadsRS
			else:
				RSCountWholeGenome = np.min([np.max(rnaSeqCountWholeGenome1), np.max(rnaSeqCountWholeGenome2)])	* uniqueReadsRSBS / uniqueReadsRS	

		

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

print('\nTS overl. known Terminators: ' + str(len(TScountsOverlappingTerminators)))
print('TS overl. genes: ' + str(len(TScountsOverlappingGenes)))
print('\n')



#######################################################################
#######################################################################
# plotting: 

npMaxTSandMeanRNAWholeGenome = np.array(maxTSandMeanRNAWholeGenome)
npTSCountsOverlappingTerminators = np.array(TScountsOverlappingTerminators)
npTScountsOverlappingGenes = np.array(TScountsOverlappingGenes)

if organism == 'B.subtilis':

# write to files with coords and strand 

	with open(outfileTSvsRS, 'w') as tsVsRs, open(outfileTScountsOverlappingTerminators, 'w') as tsOt, \
		open(outfileTScountsOverlappingGenes, 'w') as tsOg:


		for list1 in strandMaxTSandMeanRNAWholeGenome:
			tsVsRs.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n') # with coords and strand


		for list1 in strandTScountsOverlappingTerminators:
			tsOt.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n')


		for list1 in strandTScountsOverlappingGenes:
			tsOg.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n')



	npMaxTSandMeanRNATerminators = np.array(maxTSandMeanRNATerminators)

	plt.figure(1, figsize =(12,6), dpi=100)

	#at known terminators
	sub1 = plt.subplot(121) 

	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(1, 100000)
	plt.ylim(1, 100000)
	plt.grid(True)
	# adding +1 to all values --> log not - infinity
	plt.scatter(npMaxTSandMeanRNATerminators[:,1]+1, npMaxTSandMeanRNATerminators[:,0]+1, s = 15, color = 'red')
	plt.title("Term-Seq vs. RNA-Seq at known terminators \n(avg of " + str(numberOfNucsToAvg) + " nucleotide intervals)", fontsize=14)
	plt.ylabel("Max. Term-Seq count", fontsize=14)
	plt.xlabel("Avg. RNA-Seq count", fontsize=14)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	sub1.set_axisbelow(True)





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
	plt.title("Term-Seq vs. RNA-Seq \n(avg of " + str(numberOfNucsToAvg) + " nucleotide intervals)", fontsize=14)
	plt.xlabel("Avg. RNA-Seq count", fontsize=14)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	sub2.set_axisbelow(True)


	plt.savefig(outfile, dpi=300)


	plt.close()


#######################################################################



# if organism == 'E.faecalis' or organism == 'L.monocytogenes' or organism == 'S.pneumoniae':
else:

	with open(outfileTSvsRS, 'w') as tsVsRs, \
		open(outfileTScountsOverlappingGenes, 'w') as tsOg:


		for list1 in strandMaxTSandMeanRNAWholeGenome:
			tsVsRs.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n') # with coords and strand


		for list1 in strandTScountsOverlappingGenes:
			tsOg.write(str(list1[0]) + '\t' + str(list1[1]) + '\t' + str(list1[2]) + '\t' + str(list1[3]) + '\n')



	with open(RNIEoutfileTSvsRS, 'w') as RNIEtsVsRs, \
		open(RNIEoutfileTScountsOverlappingGenes, 'w') as RNIEtsOg:


		for list1 in strandMaxTSandMeanRNAWholeGenome:
			if list1[0]-l2 > 0 and list1[0]+l2 <= lengthGenome:
				if list1[3] == '-':
					RNIEtsVsRs.write(chrom + '\t' + str(list1[0] -l2) + '\t' + str(list1[0]+l1) + '\t' \
						+ str(list1[0]) + '_'+ str(list1[1]) + '_'+ str(list1[2]) + '_' + str(list1[3]) +'\n') 
				if list1[3] == '+':
					RNIEtsVsRs.write(chrom + '\t' + str(list1[0] - l1) + '\t' + str(list1[0]+l2) + '\t' \
						+ str(list1[0]) + '_' + str(list1[1]) + '_'+ str(list1[2]) + '_' +str(list1[3]) +'\n') 
				if list1[3] == 'nostrand':
					RNIEtsVsRs.write(chrom + '\t' + str(list1[0] -l2) + '\t' + str(list1[0]+l1) + '\t' \
						+ str(list1[0]) + '_' + str(list1[1]) + '_'+ str(list1[2])+ '_' +str(list1[3]) +'\n') 
					RNIEtsVsRs.write(chrom + '\t' + str(list1[0] - l1) + '\t' + str(list1[0]+l2) + '\t' \
						+ str(list1[0]) + '_' + str(list1[1]) + '_'+ str(list1[2])+ '_'  +str(list1[3]) +'\n') 

		for list1 in strandTScountsOverlappingGenes:
			if list1[0]-l2 > 0 and list1[0]+l2 <= lengthGenome:
				if list1[3] == '-':
					RNIEtsOg.write(chrom + '\t' + str(list1[0] -l2) + '\t' + str(list1[0]+l1) + '\t' \
						+ str(list1[0]) + '_' + str(list1[1]) + '_'+ str(list1[2])+ '_'  +str(list1[3]) +'\n')
				if list1[3] == '+':
					RNIEtsOg.write(chrom + '\t' + str(list1[0] - l1) + '\t' + str(list1[0]+l2) + '\t' \
						+ str(list1[0]) + '_' + str(list1[1]) + '_'+ str(list1[2])+ '_'  +str(list1[3]) +'\n')
				if list1[3] == 'nostrand':
					RNIEtsOg.write(chrom + '\t' + str(list1[0] -l2) + '\t' + str(list1[0]+l1) + '\t' \
						+ str(list1[0]) + '_' + str(list1[1]) + '_'+ str(list1[2])+ '_'  +str(list1[3]) +'\n')
					RNIEtsOg.write(chrom + '\t' + str(list1[0] - l1) + '\t' + str(list1[0]+l2) + '\t' \
						+ str(list1[0]) + '_' + str(list1[1]) + '_'+ str(list1[2])+ '_' +str(list1[3]) +'\n')




	plt.figure(1, figsize =(6,6), dpi=100)
	ax = plt.axes()

	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(1, 100000)
	plt.ylim(1, 100000)

	plt.scatter(npMaxTSandMeanRNAWholeGenome[:,2]+1, npMaxTSandMeanRNAWholeGenome[:,1]+1, s = 15, label = "non-overlapping")
	plt.scatter(npTScountsOverlappingGenes[:,2]+1, npTScountsOverlappingGenes[:,1]+1, s = 15, color = 'plum', label = "overlapping genes")


	plt.legend(prop={'size': 14})
	plt.grid(True)
	plt.title("Term-Seq vs. RNA-Seq \n(avg of " + str(numberOfNucsToAvg) + " nucleotide intervals)", fontsize=14)
	plt.xlabel("Avg. RNA-Seq count", fontsize=14)
	plt.xticks(fontsize=14)
	plt.yticks(fontsize=14)

	ax.set_axisbelow(True)


	plt.savefig(outfile, dpi=300)


	plt.close()
