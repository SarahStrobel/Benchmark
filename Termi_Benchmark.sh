# !/bin/bash

# set for Python2
PYTHON=python
SCRIPTS=Scripts-python2

# comment the following line out to disable the more time-consuming tests with 100 negatives per positive.  un-comment it to enable 100 negatives per positive
#ENABLE_100_NEG_PER_POS=1

# die if there's an error
set -e
set -o pipefail
set -x

# test if cmsearch is available, and version 1.0 -- otherwise RNie will call the wrong function
echo testing that cmsearch is version 1.0
cmsearch -h > /dev/null # check if cmsearch is available at all.  if not, you should install infernal 1.0
if cmsearch -h | grep -q "INFERNAL 1.1" ; then
    echo The cmsearch command from Infernal version 1.1 is in the PATH.  You need to make sure that version 1.0 is in the PATH, so that the RNie program will work correctly
    exit 1
fi
if cmsearch -h | grep -q "INFERNAL 1.0" ; then
    echo good, got version 1.0
else
    echo Found cmsearch command, but it is not from Infernal version 1.0.  Is this some super-old version of it?  RNie requires version 1.0
fi

# test if other required software is available
echo testing that required software is installed
rnamotif -v > /dev/null
$PYTHON -c "import argparse, collections, csv, glob, itertools, linecache, math, matplotlib, matplotlib.pyplot, numpy, os, os.path, pandas, re, subprocess, sys"  > /dev/null
$PYTHON -c "from __future__ import absolute_import, division, print_function" > /dev/null
$PYTHON -c "from bisect import bisect_left" > /dev/null
$PYTHON -c "from collections import Counter" > /dev/null
$PYTHON -c "from matplotlib.colors import ListedColormap" > /dev/null
$PYTHON -c "from operator import itemgetter" > /dev/null
$PYTHON -c "from tabulate import tabulate" > /dev/null
python3 -c "import itertools, numpy, os, pandas, subprocess, sys"
echo GOOD, looks like all required software is available


pathToParentDirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

mkdir -p $pathToParentDirectory/"Benchmark"

brev=('BS' 'EF_chrom' 'LM' 'SP' )

mkdir -p $pathToParentDirectory/Benchmark/RNIE/
mkdir -p $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive
mkdir -p $pathToParentDirectory/Benchmark/RNAmotif/
mkdir -p $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive
mkdir -p $pathToParentDirectory/Benchmark/iTerm_PseKNC/
mkdir -p $pathToParentDirectory/Benchmark/iTerm_PseKNC/1NegativePerPositive

if [[ $ENABLE_100_NEG_PER_POS == 1 ]] ; then
mkdir -p $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive
mkdir -p $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive
fi


for file in $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/*.fasta
do
	echo $file
	if /usr/bin/time $pathToParentDirectory/Termi/RNIE/rnie.pl \
		--gene \
		-f $file \
	 	-md $pathToParentDirectory/Termi/RNIE/models\
		-th 0 \
	 	-p $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/1_predicted_"$(basename "$file" .fasta)"_th0; then
	 		echo 'worked'
	else
		exit 1
	fi
	if /usr/bin/time rnamotif \
		-descr $pathToParentDirectory/Termi/RNAmotif/terminator-lesnik.desc \
		$file | \
		rmfmt -l | \
		rmprune > $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/1_predicted_"$(basename "$file" .fasta)".terminator-lesnik.out; then
			echo 'worked'
	else 
		exit 1
	fi
done

if [[ $ENABLE_100_NEG_PER_POS == 1 ]] ; then

for file in $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/*.fasta
do
	echo $file
	if /usr/bin/time $pathToParentDirectory/Termi/RNIE/rnie.pl \
		--gene \
		-f $file \
	 	-md $pathToParentDirectory/Termi/RNIE/models\
		-th 0 \
	 	-p $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive/100_predicted_"$(basename "$file" .fasta)"_th0; then
	 		echo 'worked'
	else
		exit 1
	fi
	if /usr/bin/time rnamotif \
		-descr $pathToParentDirectory/Termi/RNAmotif/terminator-lesnik.desc \
		$file | \
		rmfmt -l | \
		rmprune > $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive/100_predicted_"$(basename "$file" .fasta)".terminator-lesnik.out; then
			echo 'worked'
	else 
		exit 1
	fi
done
fi # 100 neg per pos

#cd $pathToParentDirectory/Termi/iTerm-PseKNC_modified/

for file in $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/*.csv
do
	echo $file
	if /usr/bin/time python3 $pathToParentDirectory/Termi/iTerm-PseKNC_modified/iTerm-PseKNC_modified.py \
		$file \
		$pathToParentDirectory/Benchmark/iTerm_PseKNC/1NegativePerPositive/1_predicted_"$(basename "$file" .csv)"; then
			echo 'worked'
	else
		exit 1
	fi
done

#####################################################################
#####################################################################



for file in $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/*.gff
do
	echo $file
	perl $pathToParentDirectory/Termi/RNIE/rnie2gff.pl -r $file
done


for file in $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/*.out
do
	perl $pathToParentDirectory/Termi/RNIE/terminator-lesnik2gff.pl -t $file 
done

fastaFiles=()

for file in $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/*_shuffled.fasta
do
	fastaFiles=("${fastaFiles[@]}" "$file")
done

resultFiles=()

for file in $pathToParentDirectory/Benchmark/iTerm_PseKNC/1NegativePerPositive/*result.txt
do
	resultFiles=("${resultFiles[@]}" "$file")
done

for ((i=0;i<${#resultFiles[@]};++i))
do
	echo ${resultFiles[i]}
	echo ${fastaFiles[i]}
	$PYTHON $pathToParentDirectory/$SCRIPTS/iTerm2gff.py \
	-iterm ${resultFiles[i]} \
	-i ${fastaFiles[i]} \
	-o $pathToParentDirectory/Benchmark/iTerm_PseKNC/iTerm_results_from_wl/"$(basename "${resultFiles[i]}" .txt)"
done

if [[ $ENABLE_100_NEG_PER_POS == 1 ]] ; then

for file in $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive/*.gff
do
	echo $file
	perl $pathToParentDirectory/Termi/RNIE/rnie2gff.pl -r $file
done


for file in $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive/*.out
do
	perl $pathToParentDirectory/Termi/RNIE/terminator-lesnik2gff.pl -t $file 
done
fi # 100 neg per pos


#######################################################################
#######################################################################

mkdir -p $pathToParentDirectory/"Benchmark/ROC"
mkdir -p $pathToParentDirectory/"Benchmark/ROC/1NegativePerPositive"

if [[ $ENABLE_100_NEG_PER_POS == 1 ]] ; then

mkdir -p $pathToParentDirectory/"Benchmark/ROC/100NegativesPerPositive"
fi

for ((i=0;i<${#brev[@]};++i))
do
	a=$(cat $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)
	b=$(cat $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)

	$PYTHON $pathToParentDirectory/$SCRIPTS/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/iTerm_PseKNC/iTerm_results_from_wl/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled_result.gff \
	-false $pathToParentDirectory/Benchmark/iTerm_PseKNC/iTerm_results_from_wl/1_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled_result.gff \
	-o $pathToParentDirectory/Benchmark/ROC/1NegativePerPositive/iTerm_${brev[i]}_shuffled_ \
	-nucs $((a+b))

	$PYTHON $pathToParentDirectory/$SCRIPTS/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled_th0-geneMode-rnie.bits.gff \
	-false $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/1_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled_th0-geneMode-rnie.bits.gff  \
	-o $pathToParentDirectory/Benchmark/ROC/1NegativePerPositive/RNIEth0_${brev[i]}_shuffled_ \
	-nucs $((a+b))

	$PYTHON $pathToParentDirectory/$SCRIPTS/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled.terminator-lesnik.out.dG_score.gff\
	-false $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/1_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled.terminator-lesnik.out.dG_score.gff \
	-o $pathToParentDirectory/Benchmark/ROC/1NegativePerPositive/RNAmotif_${brev[i]}_shuffled_ \
	-nucs $((a+b))
done

if [[ $ENABLE_100_NEG_PER_POS == 1 ]] ; then

for ((i=0;i<${#brev[@]};++i))
do
	a=$(cat $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)
	b=$(cat $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)

	$PYTHON $pathToParentDirectory/$SCRIPTS/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled_th0-geneMode-rnie.bits.gff \
	-false $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive/100_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled_th0-geneMode-rnie.bits.gff \
	-o $pathToParentDirectory/Benchmark/ROC/100NegativesPerPositive/RNIEth0_${brev[i]}_shuffled_ \
	-nucs $((a+b))

	$PYTHON $pathToParentDirectory/$SCRIPTS/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled.terminator-lesnik.out.dG_score.gff\
	-false $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive/100_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled.terminator-lesnik.out.dG_score.gff \
	-o $pathToParentDirectory/Benchmark/ROC/100NegativesPerPositive/RNAmotif_${brev[i]}_shuffled_ \
	-nucs $((a+b))
done
fi # 100 neg per pos

#######################################################################
#######################################################################

echo done all tests
