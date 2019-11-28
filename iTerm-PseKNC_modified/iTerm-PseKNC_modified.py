# -*- coding: utf-8 -*-
"""
Created on 2019-7-16

@author: feng
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys
import subprocess
import os
import itertools
import pandas as pd
import numpy as np

def SeqToSubSeq(seq, seq_id, slideFile):
    i = 1
    loopNum = len(seq) - 80
    
    if (loopNum <= 0):
        slideFile.write(">" + '\t' + str(seq_id) + '\t' + str(i) + '\n')
        slideFile.write(seq) 
    while (loopNum >= 1):
        slideFile.write(">" + '\t' + str(seq_id) + '\t' + str(i) + '\n')
        slideFile.write(seq[i-1:i+79] + '\n')
        i += 1
        loopNum -= 1;

def ParseSeq(seqFileName, slideFileName):
    seqFile = open(seqFileName, "r")
    slideFile = open(slideFileName, "w")
    
    seqs = seqFile.readlines()
    seqNum = len(seqs)
    for i in range(seqNum):
        if (seqs[i][0] == ">"):
            SeqToSubSeq(seqs[i+1], (i+2) // 2, slideFile)
        i += 2
    
    seqFile.close()
    slideFile.close()
 
def Bio_list(file):
    f = open(file, "r")
    lines = f.readlines()
    
    bio_list = []
    for line in lines:
        l = line.strip()
        bio_list.append(int(l))
        
    return bio_list

def ObtainKnucleotides(k):
    bases = ['A', 'T', 'C', 'G']
    k_bases = []
    k_nucleotides = []
    indexs = [''.join(x) for x in itertools.product('0123', repeat=k)]  #生成0123的排列组合
    
    for i in range(k):
        k_bases.append(bases)
    
    for index in indexs:
        k_indexs = list(index)
        m = ''
        for k_index in k_indexs:
            m = m + k_bases[k_indexs.index(k_index)][int(k_index)]
        k_nucleotides.append(m)
    
    return k_nucleotides  

def ObtainSequences(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    sequences = []
    
    for line in lines:
        if line[0] != '>':
            each_line = line.strip()
            sequences.append(each_line.upper())
    
    return sequences

def obtainNucleotidesPhysicoChemicalDict(filename):
    phychemdict = dict()
    count_line = 0
    f = open(filename)
    for eachline in f:
        count_line += 1
        temp = eachline.strip().split("\t")
        if count_line == 1:
            nucleotides = temp[0:]
        else:
            phychemdict[temp[0]] = dict()
            count_num = 0
            for each in temp[1:]:
                phychemdict[temp[0]][nucleotides[count_num]] = float(each)
                count_num += 1
    f.close()

    return phychemdict

def calculateOccurenceFrequencyOfOlineucletide(kTuples, k_nucleotides):
    occurfrequency = []
    tupleLen = len(kTuples)
    for each in k_nucleotides:
        occurfrequency.append(kTuples.count(each)/tupleLen)

    return occurfrequency

def calculateFeatureValueByCorrFactsDictAndOccurfrequency(corrFactorsDict, occurfrequency, bio_list, lamdas, nucleos, phyChemName):
    select_fy = []
    for i in bio_list:
        select_fy.append(occurfrequency[i]) 
    
    for phyChemKey in phyChemName:
        for i in range(lamdas):
            select_fy.append(corrFactorsDict[phyChemKey][i+1])
    
    return select_fy

def fn_physic(sequence, k, k_nucleotides, nucleStandDict, bio_list, lamdas):
    corrFactorsDict = dict()
    seqLen = len(sequence)

    kTuples = []
    for i in range(seqLen - k + 1):
        kTuples.append(sequence[i:i+k])
    occurfrequency = calculateOccurenceFrequencyOfOlineucletide(kTuples, k_nucleotides)

    kTuples = []
    for i in range(seqLen - 2 + 1):
        kTuples.append(sequence[i:i+2])

    phyChemName = list(nucleStandDict.keys())
    for eachName in phyChemName:
        corrFactorsDict[eachName] = dict()
        for lamba in range(1, (lamdas + 1)):
            temp = []
            for kTuplesIndex in range(len(kTuples)-lamba):
                preKTuple = kTuples[kTuplesIndex]
                backKTuple = kTuples[kTuplesIndex+lamba]
                tempNumber = nucleStandDict[eachName][preKTuple]*nucleStandDict[eachName][backKTuple]
                temp.append(tempNumber)
            corrFactorsDict[eachName][lamba] = np.mean(temp)

    pse = calculateFeatureValueByCorrFactsDictAndOccurfrequency(corrFactorsDict, occurfrequency, bio_list, lamdas, k_nucleotides, phyChemName)    

    return pse

def iterm_pseknc_process(bioFile, inputFile, outputFile):
    bio_list = Bio_list(bioFile)
    k_nucleotides = ObtainKnucleotides(5)
    sequences = ObtainSequences(inputFile)
    nucleStandDict = obtainNucleotidesPhysicoChemicalDict(sys.argv[3] + "6_standard.txt")
    outputFile = open(outputFile, "w")
    for sequence in sequences:
        pse = fn_physic(sequence, 5, k_nucleotides, nucleStandDict, bio_list, 5)
        
        outputFile.write("-1")
        for i in range(len(pse)):	
            outputFile.write("\t%d:%.6f" % (i+1,pse[i]))
        outputFile.write("\n") 

def read_result(slideFileName, PredictFile, outputFileName): 
    slideFile = open(slideFileName, 'r')
    predictFile = open(PredictFile, 'r')
    outputFile = open(outputFileName, 'w')
    
    slidlines = slideFile.readlines()
    resultlines = predictFile.readlines()[1:]
    slidnum = len(slidlines)
    resulnum = len(resultlines)
    for i in range(resulnum):
        if (slidlines[2*i][0] == '>'):
            result = resultlines[i].strip().split(' ')
            if result[0] == '1':
                outputFile.write(slidlines[2*i][0:-1] + " is a terminator, the possible is " + str(result[1]) + '\n')
            if result[0] == '-1':
                outputFile.write(slidlines[2*i][0:-1] + " is a non-terminator, the possible is " + str(result[2]) + '\n')  
            outputFile.write(slidlines[2*i+1])
   
    slideFile.close()
    predictFile.close()
    outputFile.close()
    
  
InputFile = sys.argv[1]
slideFile = sys.argv[2] + "_SlideFile.txt"
bioFile = sys.argv[3] + "Bio_index.txt"
ToScaleFile = sys.argv[2] + "_ToScaleFile.txt"
ScaleFile = sys.argv[2] + "_ScaleFile.scale"
PredictFile = sys.argv[2]+ "_PredictFileName.txt"
resultFile = sys.argv[2] + "_result.txt"

if __name__ == "__main__":
    ParseSeq(InputFile, slideFile)
    iterm_pseknc_process(bioFile, slideFile, ToScaleFile)
    os.popen("svm-scale -r iterm.rule ToScaleFile.txt > ScaleFile.scale")
    os.popen("svm-predict -b 1 ScaleFile.scale iterm.model PredictFileName.txt")
    read_result(slideFile, PredictFile, resultFile)
    os.remove(slideFile)
    os.remove(ToScaleFile)
    os.remove(ScaleFile)
    os.remove(PredictFile)
    
