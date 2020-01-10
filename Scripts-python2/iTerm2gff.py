import argparse
import os.path
import linecache



parser = argparse.ArgumentParser(description= 'embed fasta sequences' + '\n'
                                'Usage:' + '\t' + 'iTerm2gff.py <options> -i -o')

#required files:
parser.add_argument('-iterm', dest='iTermOutput', help='output file from iTerm-PseKNC', required=True)
parser.add_argument('-i', dest='iTermInput', help='input file for iTerm-PseKNC in FASTA format', required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)

args = parser.parse_args()

iTermOutput = args.iTermOutput
iTermInput = args.iTermInput
outpath = args.outpath


outfile = outpath + '.gff'

ids = []

#######################################################################
#######################################################################

# id	iTerm-PseKNC	terminator	start	end		score	strand	frame (bei rnie2gff einfach -)	attribute (additional info)

with open(iTermOutput, 'r') as iTermOut, open(iTermInput, 'r')  as iTermIn, open(outfile, 'w') as out:

	for line in iTermIn:
		if '>' in line:
			ids.append(line[1:].strip())



	for line  in iTermOut:

		if 'is a terminator' in line:
			IDiTerm = int(line.split()[1])
			ID = ids[IDiTerm-1]
			start = int(line.split()[2])		
			score = float(line.split()[-1])
			sequence = iTermOut.next()	
			end = start + (len(sequence)-1)

			out.write(str(ID) + '\t' + 'iTerm-PseKNC' + '\t' + 'terminator' + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\n')

		if 'non-terminator' in line:
			IDiTerm = int(line.split()[1])
			ID = ids[IDiTerm-1]
			start = int(line.split()[2])		
			score = 1-float(line.split()[-1])
			sequence = iTermOut.next()	
			end = start + (len(sequence)-1)
			
			out.write(str(ID) + '\t' + 'iTerm-PseKNC' + '\t' + 'non-terminator' + '\t' + str(start) + '\t' + str(end) + '\t' + str(score) + '\t' + '.' + '\t' + '.' + '\t' + '.' + '\n')

