# counting the 5' positions of each term-seq read and creating a bedgraph file
import glob
import numpy as np
import argparse
import re
from collections import Counter
from collections import OrderedDict

#######################################################################
#######################################################################
class OrderedCounter(Counter, OrderedDict):
    'Counter that remembers the order elements are first encountered'

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, OrderedDict(self))

    def __reduce__(self):
        return self.__class__, (OrderedDict(self),)



#######################################################################
#######################################################################

parser = argparse.ArgumentParser(description= 'Count 5\' end nuc positions' + '\n'
								'Usage:' + '\t' + '5primePosfromSam_2_Bedgraph.py <options> -rnaSeq -termSeq -o')

#required files:
parser.add_argument('-rnaSeq', dest='rnaSeqFile', help='input of aligned RNA-Seq coverage file in bedgraph format', required=True)
parser.add_argument('-termSeq', dest='termSeqPath', help='path to input of aligned Term-Seq alignement file in sam format', required=True)
parser.add_argument('-o', dest='outfile', help='output path and filename prefix', required=True)


args = parser.parse_args()

rnaSeqFile = args.rnaSeqFile
outfile = args.outfile
termSeqPath = args.termSeqPath



#Bacillus subtilis
chrom = "gi|255767013|ref|NC_000964.3|"
lengthGenome = 4215606


#######################################################################
#######################################################################

for samFile in glob.glob(termSeqPath+'*.sam'):

	experiment = samFile.split('/')[-1].split('.')[0]
	print experiment

	fivePrimeEndsWithStrands = []
	coordDictWithStrands = {}


	with open(samFile, 'r') as f1:
		with open(outfile + experiment + '.bedgraph', 'w') as ofBEDGRAPH:
		
			for line in f1:
				if '@' not in line:			
					QNAME = line.split()[0]
					FLAG = int(line.split()[1])
					POS = int(line.split()[3])
					CIGAR = line.split()[5]
					SEQ = line.split()[9]

					end =  POS + len(SEQ)-1
					strand = ''

					matches = 0

					# http://broadinstitute.github.io/picard/explain-flags.html
					# 1 (0X1): read paired, 4 (0x4): read unmapped, 8 (0X8): mate unmapped, 
					# 16 (0x10): read reverse strand, 32 (0x20): mate reverse strand,
					# 64 (0x40): first in pair, 128 (0x80): second in pair,
					# 256 (0x100): not primary alignment, 272 (0x100, 0x10): not primary alignment reverse strand,
					# 516 (0x200): read fails platform/vendor quality check, 1024 (0x400): read is pcr or optical duplicate,
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


			coordDictWithStrands = Counter(fivePrimeEndsWithStrands)

			newCoordDictWithStrands = {}

			for key, value in coordDictWithStrands.items():
				newCoordDictWithStrands[key[0]] = [value, key[1]]

#######################################################################
#######################################################################

			for start in range(1,lengthGenome+1):
				if start in newCoordDictWithStrands:
					ofBEDGRAPH.write(str(chrom) + '\t' + str(start) + '\t' + str(start+1) + '\t' +str(newCoordDictWithStrands[start][0]) + '\t' + str(newCoordDictWithStrands[start][1]) + '\n')
				else:
					ofBEDGRAPH.write(str(chrom) + '\t' + str(start) + '\t' + str(start+1) + '\t' + str(0) + '\t' + 'nostrand' + '\n')


		ofBEDGRAPH.close()
		f1.close()



