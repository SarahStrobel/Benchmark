import math
import csv
from operator import itemgetter
import argparse
import os.path
import itertools

#######################################################################
#######################################################################
# methods for checking parsed file types

def checkFastaFormat(v):
    b = os.path.splitext(v)[1][1:].lower()
    if not (b == 'fa' or b == 'fasta'):   
        raise argparse.ArgumentTypeError('fasta format file type expected')
    else:
        return v

def checkGffFormat(v):
    b = os.path.splitext(v)[1][1:].lower()
    if b != 'gff':
        raise argparse.ArgumentTypeError('gff format file type expected')
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
# #overlap function

def overlap(coords1,coords2):   
    coords1.sort()
    coords2.sort()

    x3 = max(coords1[0],coords2[0]) #start1, start2
    y3 = min(coords1[1],coords2[1]) #end1, end2

    if x3 <= y3:  
        return True
    else:
        return False


#######################################################################
# #functions for sens, ppv, mcc

def calcSens(tp, fn):
    denom = float(tp + fn)
    sens = 0.00
    if (denom > 0):
        sens = tp / denom
        if sens > 1.0:
            sens = 1.0
    
    return sens



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

parser = argparse.ArgumentParser(description= 'compute ROC values' + '\n'
                                'Usage:' + '\t' + 'computeROC.py <options> -pos -neg -true -false -o')

#required files:
parser.add_argument('-pos', dest='positivesFile', help='input predicted Terminators', type=checkFastaFormat, required=True)
parser.add_argument('-neg', dest='negativesFile', help='input predicted Terminators', type=checkFastaFormat, required=True)
parser.add_argument('-true', dest='trueFile', help='input outfiles from prediction tool', type=checkGffFormat, required=True)
parser.add_argument('-false', dest='falseFile', help='input outfiles from prediction tool', type=checkGffFormat, required=True)
parser.add_argument('-o', dest='outpath', help='output path and filename prefix', required=True)
parser.add_argument('-nucs', dest='totalNucs', help='total number of nucleotides in positive and negative FASTA files', required=True)
parser.add_argument('-l', dest='lengthTerminator', help='length of artificial terminator, default:100', type=checkInt, nargs='?', default=100)
parser.add_argument('-n', dest='lengthOfRegion', help='length of region spanning the termination signal', type=checkInt, nargs='?', default=40)

args = parser.parse_args()

positivesFile = args.positivesFile
negativesFile = args.negativesFile
trueFile = args.trueFile
falseFile = args.falseFile
outpath = args.outpath
totalNucs = int(args.totalNucs)
lengthOfRegion = int(args.lengthOfRegion)
lengthTerminator = int(args.lengthTerminator)

l1 = lengthTerminator * 0.16666666666666664
l2 = lengthTerminator - l1

l1 =  int(math.ceil(l1))
l2 = int(math.ceil(l2))


l3 = lengthOfRegion * 0.25
l4 = lengthOfRegion - l3

l3 =  int(math.ceil(l3))
l4 = int(math.ceil(l4))


print totalNucs

outfile = outpath + 'terminator.out.dG_score-accuracy.dat'


print outfile


#######################################################################
#######################################################################

numPositives = 0
numNegatives = 0

numOfPosPred = 0
numOfNegPred = 0

# dictionaries that save every true/false terminator from predicted positives and predicted negatives FASTA files (key = id, value = 'T' bzw 'F', score)
trueDict = {}
falseDict = {}

tntpList = []

#placeholder of 0, if input smaller than 0 change
minScore = 0 
numberOfUniquePos=0
numberOfUniqueNeg=0

#read input predicted positives fasta file and put into trueDict
with open(positivesFile, 'r') as f1, open(negativesFile, 'r') as f2, \
    open(trueFile, 'r') as f3, open(falseFile, 'r') as f4:

    for line in f1:  
        words = line.split()

        if '>' in words[0]:
            identifier = words[0]
            strand = words[0].split('_')[-1]

            # embedded
            start = 501 + l2 - l4
            end = (501 + lengthTerminator) - l1 + l3 

            coords = [start, end]
        
            trueElement = ['T', minScore, coords, identifier[1:]]

            #check if identifier is only used once, put in trueDict
            key = identifier[1:]
            if key in trueDict:
                raise Exception("multiple entries with same ID in predicted positives fasta file")
            else:
                trueDict[key] = trueElement

            numPositives += 1

    

    # #read input shuffled negatives fasta file and put into falseDict

    for line in f2:
        words = line.split()
        
        if '>' in words[0]:
            identifier = words[0]
            strand = words[0].split('_')[-1]

            # embedded
            start = 501 + l2 - l4
            end = (501 + lengthTerminator) - l1 + l3 

            coords = [start, end]

            falseElement = ['F', minScore, coords, identifier[1:]]

            key = identifier[1:] #slice to omit first character ('>')
            if key in trueDict:
                raise Exception("multiple entries with same ID in predicted negatives fasta file")
            else:
                falseDict[key] = falseElement

            numNegatives += 1



    # #read output of terminator prediction software for predicted positives

    for line in f3:
        identifier = line.split()[0]
        start = int(line.split()[3])
        end = int(line.split()[4])
        coords = [start, end]
        score = float(line.split()[5])
        numOfPosPred += 1

        if score < minScore:
            raise Exception("score in rnieTrue < minScore; change minScore")


        #check if overlap
        if overlap(trueDict[identifier][2], coords):
        #if score in dict > new score --> ignore, else replace with higher score
            if trueDict[identifier][1] < score:
                trueDict[identifier][1] = score


                tntpList.append(trueDict[identifier])

    tntpList.sort()
    tntpList=list(tntpList for tntpList,_ in itertools.groupby(tntpList))
    # print tntpList
    numberOfUniquePos = len(tntpList)
    print 'number of unique positives: ' + str(numberOfUniquePos)


    # #read output of terminator prediction software for predicted negatives
    for line in f4:
        identifier = line.split()[0]
        start = int(line.split()[3])
        end = int(line.split()[4])
        coords = [start, end]
        score = float(line.split()[5])
        numOfNegPred += 1


        if score < minScore:
            raise Exception("score in rnieFalse < minScore; change minScore")



#removed overlap because it was too different from Paul's solution

        # if overlap(falseDict[identifier][2], coords):
        #     if falseDict[identifier][1] < score:
        #         falseDict[identifier][1] = score

        #         tntpList.append(falseDict[identifier])

        if falseDict[identifier][1] < score:
            falseDict[identifier][1] = score

            tntpList.append(falseDict[identifier])

    tntpList.sort()
    tntpList=list(tntpList for tntpList,_ in itertools.groupby(tntpList))
    numberOfUniqueNeg = len(tntpList) - numberOfUniquePos
    print 'number of unique negatives: ' + str(numberOfUniqueNeg)


print "number of positives: " + str(numPositives)
print "number of negatives: " + str(numNegatives)
print "positives predicted by tool: " + str(numOfPosPred)
print "negatives predicted by tool: " + str(numOfNegPred)

# print tntpList

tntpList.sort(key=itemgetter(1),reverse=True)


# getting tp,fn,fp,tn for each element to compute ROC curve
# starting at +infinite

tp=0
fp=0
tn=numNegatives
fn=numPositives


lastSens=0
lastPPV=0
lastMCC=0
lastScore=minScore


with open(outfile,"w") as of:
    
    of.write(str("sens" + "\t" + "ppv"+ "\t" + "mcc" + "\t" + "fpr/KB" + "\t" + "score")+"\n")

    for score in tntpList:
        if "T" in score[0]:
            tp+=1 
            fn-=1 

        else:
            fp+=1 
            tn-=1 


        sens = calcSens(tp,fn)
        ppv = calcPpv(tp,fp)
        mcc = calcMcc(tp,tn,fp,fn)
        fpr = (fp/float(totalNucs))*1000
     

        if ((abs(lastSens-sens) > 0.001 or abs(lastPPV-ppv)>0.001 or abs(lastMCC-mcc)>0.001) and lastScore != score[1]):
            toWrite = str("{:.6f}".format(sens)) +"\t"+ str("{:.6f}".format(ppv)) +"\t"+ str("{:.6f}".format(mcc)) + "\t" + "{:.6f}".format(fpr) + "\t" + "{:.6f}".format(score[1])
            of.write(toWrite + "\n")

            lastSens=sens
            lastPPV=ppv
            lastMCC=mcc
            lastScore=score[1]

