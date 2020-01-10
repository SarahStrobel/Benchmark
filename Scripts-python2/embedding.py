import argparse
import os.path


#######################################################################
#######################################################################
# methods for checking parsed file types

def checkFastaFormat(v):
    b = os.path.splitext(v)[1][1:].lower()
    if not (b == 'fa' or b == 'fasta'):   
        raise argparse.ArgumentTypeError('fasta format file type expected')
    else:
        return v

#######################################################################
#######################################################################

parser = argparse.ArgumentParser(description= 'embed fasta sequences' + '\n'
                                'Usage:' + '\t' + 'embedding.py <options> -term -neg -front -back -frontN -backN -o')

#required files:
parser.add_argument('-term', dest='predictedTerminatorFile', help='sequences of predicted terminators', type=checkFastaFormat, required=True)
parser.add_argument('-neg', dest='predictedNegativeFile', help='shuffled sequences negatives', type=checkFastaFormat, required=True)
parser.add_argument('-front', dest='frontFile', help='shuffled sequences in front of predicted terminators', type=checkFastaFormat, required=True)
parser.add_argument('-back', dest='backFile', help='shuffled sequences behind predicted terminators', type=checkFastaFormat,required=True)
parser.add_argument('-frontN', dest='nativeFrontFile', help='native sequences in front of predicted terminators', type=checkFastaFormat, required=True)
parser.add_argument('-backN', dest='nativeBackFile', help='native sequences behind predicted terminators', type=checkFastaFormat,required=True)

parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)

args = parser.parse_args()

predictedTerminatorFile = args.predictedTerminatorFile
predictedNegativeFile = args.predictedNegativeFile
frontFile = args.frontFile
backFile = args.backFile
nativeFrontFile = args.nativeFrontFile
nativeBackFile = args.nativeBackFile
outpath = args.outpath


outfile1 = outpath + 'embedded_predictedTerminators_native.fasta'
outfile2 = outpath + 'embedded_predictedTerminators_shuffled.fasta'
outfile3 = outpath + 'embedded_shuffledNegatives_shuffled.fasta'





with open(predictedTerminatorFile, 'r') as term, \
	open(frontFile, 'r') as front, open(backFile, 'r') as back, \
	open(nativeFrontFile, 'r') as frontN, open(nativeBackFile, 'r') as backN, \
	open(outfile1, 'w') as out1, open(outfile2, 'w') as out2:


	for line1, line2, line3, line4, line5 in zip(frontN,term,backN,front,back):
		#native front and back + native predicted terminators
		if '>' in line1:
			out1.write(line1)
		else:
			out1.write(str(line1.rstrip())+str(line2.rstrip())+str(line3.rstrip())+"\n")
		#shuffled front and back + native predicted terminators
		if '>' in line2:
			out2.write(line2)
		else:
			out2.write(str(line4.rstrip())+str(line2.rstrip())+str(line5.rstrip())+"\n")


with open(predictedNegativeFile, 'r') as neg,\
	open(frontFile, 'r') as front2, open(backFile, 'r') as back2, \
	open(outfile3, 'w') as out3:

	for line6, line7, line8 in zip(front2,neg,back2):
		#shuffled front and back + shuffled negatives
		if '>' in line7:
			
			out3.write(line7)
		else:
			out3.write(str(line6.rstrip())+str(line7.rstrip())+str(line8.rstrip())+"\n")









