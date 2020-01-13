# counting the 5' positions of each term-seq read and creating a bedgraph file
import glob
import numpy as np
import argparse
import re
from collections import Counter


#######################################################################
#######################################################################

parser = argparse.ArgumentParser(description= 'Count 5\' end nuc positions' + '\n'
								'Usage:' + '\t' + '5primePosfromSam_2_Bedgraph.py <options> -termSeq -o')

#required files:
# parser.add_argument('-rnaSeq', dest='rnaSeqFile', help='input of aligned RNA-Seq coverage file in bedgraph format', required=True)
parser.add_argument('-termSeq', dest='termSeqPath', help='path to input of aligned Term-Seq alignement file in sam format', required=True)
parser.add_argument('-o', dest='outfile', help='output path and filename prefix', required=True)


args = parser.parse_args()


outfile = args.outfile
termSeqPath = args.termSeqPath


# rnaSeqFile = args.rnaSeqFile
# rna = open(rnaSeqFile, 'r')
# lengthGenome = len(rna.readlines())


# if 'SP' in rnaSeqFile:
# 	organism = 'S.pneumoniae'
# 	chrom = "NC_003028.3"
# 	chrom2 = chrom
# 	lengthGenome = 



#######################################################################
#######################################################################
# read all samfiles in folder and count 5' end nucs


for samFile in glob.glob(termSeqPath + '*.sam'):
	organism = ''
	chrom = ''
	chrom2 = ''
	plasmid = ''
	lengthGenome = 0

	if "Bacillus" in samFile:
		organism = 'Bacillus_subtilis'
		chrom = 'NC_000964.3/1-4215606'
		chrom2 = 'gi|255767013|ref|NC_000964.3|'
		lengthGenome = 4215606
	if "Enterococcus" in samFile:
		organism = 'Enterococcus_faecalis'
		if 'chromosome' in samFile:
			chrom = "NC_004668.1"
			chrom2 = chrom
			plasmid = 'chromosome'
			lengthGenome = 3218031
		if 'plasmid1' in samFile:
			chrom = "NC_004669.1"
			chrom2 = chrom
			plasmid = 'plasmid1'
			lengthGenome = 66320
		if 'plasmid2' in samFile:
			chrom = "NC_004671.1"
			chrom2 = chrom
			plasmid = 'plasmid2'
			lengthGenome = 57660
		if 'plasmid3' in samFile:
			chrom = "NC_004670.1"
			chrom2 = chrom
			plasmid = 'plasmid3'
			lengthGenome = 17963
	if "Listeria" in samFile:
		organism = 'Listeria_monocytogenes'
		chrom = "NC_003210.1"
		chrom2 = chrom
		lengthGenome = 2944528
	if "Streptococcus" in samFile:
		organism = 'Streptocoocus_pneumoniae'
		chrom = "NC_003028.3"
		chrom2 = chrom
		lengthGenome = 2160842

	experiment = samFile.split('/')[-1].split('.')[0].split('_')[2]

	print experiment

	if organism == 'Enterococcus_faecalis':
		outfile = outfile + experiment + '_' + organism + '_' + plasmid
	else:
		outfile = outfile + experiment + '_' + organism

	print '\n' + str(organism) + ' ' + str(plasmid)
	print 'genome length: ' + str(lengthGenome)


	fivePrimeEndsWithStrands = []
	coordDictWithStrands = {}
	uniqueReads = 0


	with open(samFile, 'r') as f1:

		with open(outfile +'.bedgraph', 'w') as ofBEDGRAPH:
			

			for line in f1:
							
				QNAME = line.split()[0]

				if '@' not in QNAME:

					FLAG = int(line.split()[1])
					POS = int(line.split()[3])
					MAPQ = int(line.split()[4])
					CIGAR = line.split()[5]
					SEQ = line.split()[9]

					if MAPQ > 0:
						uniqueReads+=1

					end =  POS + len(SEQ)-1
					strand = ''
					matches = 0

					Ss = re.findall(r'(\d+)S', CIGAR)
					if not Ss:
					# http://broadinstitute.github.io/picard/explain-flags.html
					# 1 (0X1): read paired, 4 (0x4): read unmapped, 8 (0X8): mate unmapped, 
					# 16 (0x10): read reverse strand, 32 (0x20): mate reverse strand,
					# 64 (0x40): first in pair, 128 (0x80): second in pair,
					# 256 (0x100): not primary alignment, 272 (0x100, 0x10): not primary alignment reverse strand,
					# 512 (0x200): read fails platform/vendor quality check, 1024 (0x400): read is pcr or optical duplicate,
					# 2048 (0x800): supplementary alignment

					# bitwise operator: only primary unpaired reads + reverse of that
						if (FLAG & 1) != 1 and (FLAG & 2) != 2 and (FLAG & 4) != 4 and (FLAG & 8) != 8 and \
							(FLAG & 32) != 32 and (FLAG & 64) != 64 and (FLAG & 128) != 128 and (FLAG & 256) != 256 and \
							(FLAG & 512) != 512 and (FLAG & 1024) != 1024 and (FLAG & 2048) != 2048:

							Ms = re.findall(r'(\d+)M', CIGAR)
							
							for match in Ms:
								matches += int(match)
							

							if matches == 0:
								matches = len(SEQ)

							end = POS + (matches)-1

							# if positive strand
							if FLAG != 16:
								end = min([POS, end])
								strand = '+'

							# if reverse strand
							if FLAG == 16:
								strand = '-'

							fivePrimeEndsWithStrands.append((end, strand))

			print 'unique reads: ' + str(uniqueReads)
			print '\n'

			coordDictWithStrands = Counter(fivePrimeEndsWithStrands)

			newCoordDictWithStrands = {}

			for key, value in coordDictWithStrands.items():
				newCoordDictWithStrands[key[0]] = [value, key[1]]

#######################################################################
#######################################################################
		# write outfile in bedgraph format (can be viewed in a genome browser as a histogram)

			ofBEDGRAPH.write("track type=bedGraph name=\"" + str(uniqueReads) + "\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n")

			for start in range(1,lengthGenome+1):
				if start in newCoordDictWithStrands:
					ofBEDGRAPH.write(str(chrom2) + '\t' + str(start-1) + '\t' + str(start) + '\t' + str(newCoordDictWithStrands[start][0]) + '\t' + str(newCoordDictWithStrands[start][1]) + '\n')
				else:
					ofBEDGRAPH.write(str(chrom2) + '\t' + str(start-1) + '\t' + str(start) + '\t' + str(0) + '\t' + 'nostrand' + '\n')
			outfile = args.outfile

		ofBEDGRAPH.close()
		f1.close()



