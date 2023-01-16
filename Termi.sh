#!/bin/bash

# set for Python2
PYTHON=python
SCRIPTS=Scripts-python2

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

# test esl-shuffle version
echo testing that esl-shuffle is found
esl-shuffle -h > /dev/null
echo test esl-shuffle is from Infernal 1.0.2
if esl-shuffle -h | grep -q "Copyright .C. 2009"; then
    echo good
else
    echo esl-shuffle is not from Infernal version 1.0.2.  This is actually okay, but you will not get the exact same results.
    exit 1
fi

# test if other required software is available
echo testing that required software is installed
wget --help > /dev/null
awk --help > /dev/null
sort --help > /dev/null
fasterq-dump --help > /dev/null # should be in recent versions of the sra-toolkit package.  Otherwise, install the pre-compiled binaries of the sra-toolkit at the NCBI site (currently https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software )
novoindex > /dev/null
blastn -h > /dev/null # package: ncbi-blast+
samtools --help > /dev/null
bedtools > /dev/null
esl-shuffle -h > /dev/null # this command comes from the Infernal software, but is not installed by default.  In the easel/miniapps subdirectory, do 'make install'
# fastp, python3 stuff
fastp > /dev/null
$PYTHON -c "import argparse, collections, csv, glob, itertools, linecache, math, matplotlib, matplotlib.pyplot, numpy, os, os.path, pandas, re, subprocess, sys"  > /dev/null
$PYTHON -c "from __future__ import absolute_import, division, print_function" > /dev/null
$PYTHON -c "from bisect import bisect_left" > /dev/null
$PYTHON -c "from collections import Counter" > /dev/null
$PYTHON -c "from matplotlib.colors import ListedColormap" > /dev/null
$PYTHON -c "from operator import itemgetter" > /dev/null
$PYTHON -c "from tabulate import tabulate" > /dev/null
python3 -c "import itertools, numpy, os, pandas, subprocess, sys"
echo GOOD, looks like all required software is available

ls 316*.fastq.gz > /dev/null # make sure the Warrier et al files are in the current directory

pathToParentDirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


################################################################################################################################################

# download scripts, known terminators and genomes from Sarah Strobel's github to new folder 
git clone https://github.com/SarahStrobel/Benchmark.git $pathToParentDirectory"/Termi/"

printf "\n##########################################################"
printf '\n all scripts, terminators and genomes downloaded'
printf "\n##########################################################\n\n"

#################################################################################################################################################

## Actually, we put this file into our repository so that we're safe with link rot
## download known terminator file from Paul Gardners's github to new folder with wget
#wget https://raw.githubusercontent.com/ppgardne/RNIE-benchmark/master/training-test/true.fa -P $pathToParentDirectory"/Termi/knownTerminators/"


## change u to ts and make uppercase in known terminators sequences
#sed 's/U/t/g' $pathToParentDirectory"/Termi/knownTerminators/true.fa" |
#awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1 " " $2}}'  \
	> $pathToParentDirectory"/Termi/knownTerminators/true_upper_UtoT.fa"

# download B.subtilis known terminator file form the Lin group homepage
wget http://lin-group.cn/server/iTerm-PseKNC/dependent_data2.csv -O $pathToParentDirectory"/Termi/knownTerminators/dependent_data2.fa"

printf "\n##########################################################"
printf '\n\t\t Known terminators downloaded'
printf "\n##########################################################\n\n"

################################################################################################################################################

# download RNIE models and rnie2gff.pl from Paul Gardner's github to new folder
git clone https://github.com/ppgardne/RNIE.git $pathToParentDirectory"/Termi/RNIE/"
wget https://raw.githubusercontent.com/ppgardne/RNIE-benchmark/master/scripts/rnie2gff.pl -P $pathToParentDirectory"/Termi/RNIE/"
wget https://raw.githubusercontent.com/ppgardne/RNIE-benchmark/master/scripts/terminator-lesnik2gff.pl -P $pathToParentDirectory"/Termi/RNIE/"

printf "\n##########################################################"
printf '\n\t\t RNIE models downloaded'
printf "\n##########################################################\n\n"

################################################################################################################################################

# download fastq files with fasterq-dump
mkdir -p $pathToParentDirectory"/Termi/RNASeq"
mkdir -p $pathToParentDirectory"/Termi/TermSeq"

#TermseqFiles=( ERX1304415 ERX1320300 ERX1320301 ERR1248401 ERR1248402 ERR1248403 ERR1248436 ERR1248437 ERR1248438 )
TermseqFiles=( ERX1304415 ERX1320300 ERX1320301 ERR1248401 ERR1248402 ERR1248403 ERR1248436 ERR1248437 ERR1248438 SRR7160964 SRR7160965 SRR7160966 SRR7160967)
for i in "${TermseqFiles[@]}"
do
	echo $i
	fasterq-dump $i --outdir $pathToParentDirectory"/Termi/TermSeq"
done

RNAseqFiles=( ERX1320302 ERR1248404 ERR1248439 )

for i in "${RNAseqFiles[@]}"
do
	echo $i
	fasterq-dump $i --outdir $pathToParentDirectory"/Termi/RNASeq"
done

####
# here you have to ask for the RNA-Seq files from Warrier et al. somehow and save the fastq files in the $pathToParentDirectory/Termi/RNASeq/ folder :)
#####

printf 'Temporary hack for Streptococcus'
cp 316-RNASeqND0min*.fastq.gz $pathToParentDirectory/Termi/RNASeq/ # if you change the file names, also change the stuff that's hardcoded in combine_RNASeq_SP.py
gzip -d $pathToParentDirectory/Termi/RNASeq/316*.fastq.gz

printf "\n##########################################################"
printf '\n    all RNA-Seq and Term-Seq experiments downloaded'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# trim fastq files (quality control and adapter trimming)

for filename in $pathToParentDirectory/Termi/RNASeq/*.fastq
do
	fastp -i $filename -o  $pathToParentDirectory/Termi/RNASeq/"trimmed_$(basename "$filename" )"
done

for filename in $pathToParentDirectory/Termi/TermSeq/*.fastq
do
	fastp -i $filename -o $pathToParentDirectory/Termi/TermSeq/"trimmed_$(basename "$filename" )"
done

printf "\n##########################################################"
printf '\n    Adapter Sequences trimmed'
printf "\n##########################################################\n\n"

################################################################################################################################################

# create index files and align fastq files to respective genomes with novoalign

for filename in $pathToParentDirectory/Termi/Genomes/*.fasta
do
	novoindex $pathToParentDirectory/Termi/Genomes/"$(basename "$filename" .fasta).index" $filename
done

mkdir -p $pathToParentDirectory/Termi/Alignments/
mkdir -p $pathToParentDirectory/Termi/Alignments/RNASeq/
mkdir -p $pathToParentDirectory/Termi/Alignments/TermSeq/

########################
# rna seq
########################
novoalign -f $pathToParentDirectory/Termi/RNASeq/trimmed_ERX1320302.fastq \
			-d $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.index \
			-o SAM > $pathToParentDirectory/Termi/Alignments/RNASeq/trimmed_ERX1320302_Bacillus_subtilis.sam

novoalign -f $pathToParentDirectory/Termi/RNASeq/trimmed_ERR1248439.fastq \
			-d $pathToParentDirectory/Termi/Genomes/Listeria_monocytogenes.index \
			-o SAM > $pathToParentDirectory/Termi/Alignments/RNASeq/trimmed_ERR1248439_Listeria_monocytogenes.sam

for filename in $pathToParentDirectory/Termi/Genomes/Enterococcus*.index
do
	novoalign 	-f $pathToParentDirectory/Termi/RNASeq/trimmed_ERR1248404.fastq \
				-d $filename \
				-o SAM > $pathToParentDirectory/Termi/Alignments/RNASeq/trimmed_ERR1248404_"$(basename "$filename" .index).sam"
done

for filename in $pathToParentDirectory/Termi/RNASeq/trimmed_316*.fastq
do
	novoalign -f $filename \
			-d $pathToParentDirectory/Termi/Genomes/Streptococcus_pneumoniae.index \
			-o SAM > $pathToParentDirectory/Termi/Alignments/RNASeq/"$(basename "$filename" .fastq)_Streptococcus_pneumoniae.sam"
done

########################
# term seq
########################

for filename in $pathToParentDirectory/Termi/TermSeq/trimmed_ERX*
do
	novoalign -f $filename \
			-d $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.index \
			-o SAM > $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename $filename .fastq)_Bacillus_subtilis.sam"
done

for filename in $pathToParentDirectory/Termi/TermSeq/trimmed_ERR124843*
do
	novoalign -f $filename \
			-d $pathToParentDirectory/Termi/Genomes/Listeria_monocytogenes.index \
			-o SAM > $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename $filename .fastq)_Listeria_monocytogenes.sam"
done

for filename in $pathToParentDirectory/Termi/Genomes/Enterococcus*.index
do
	for filename2 in $pathToParentDirectory/Termi/TermSeq/trimmed_ERR124840*
	do
		novoalign 	-f $filename2 \
					-d $filename \
					-o SAM > $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename "$filename2" .fastq)"_$(basename "$filename" .index)".sam"
	done
done

for filename in $pathToParentDirectory/Termi/TermSeq/trimmed_SRR*
do
	echo $filename
	novoalign -f $filename \
			-d $pathToParentDirectory/Termi/Genomes/Streptococcus_pneumoniae.index \
			-o SAM > $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename "$filename" .fastq)_Streptococcus_pneumoniae.sam"
done

########################
# known terminators to B.subtilis
########################

novoalign 	-f $pathToParentDirectory/Termi/knownTerminators/true_upper_UtoT.fa \
			-d $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.index \
			-o SAM > $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.sam"


samtools view -bS $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.sam" > $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.bam"

bedtools bamtobed -i $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.bam" \
					> $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.bed"

printf "\n##########################################################"
printf '\n    all Files aligned'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# remove multimapped reads from sam files and use samtools view to convert into bam files, sort and index (viewable in igv)

for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/*
do
	echo $filename 
	awk '/^@/ || $5 > "0"' $filename > $pathToParentDirectory/Termi/Alignments/RNASeq/"filtered_$(basename "$filename")"
done


for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/filtered_*
do
	samtools view -bS $filename |
	samtools sort - -o $pathToParentDirectory/Termi/Alignments/RNASeq/"sorted_$(basename "$filename" .sam).bam"  > $pathToParentDirectory/Termi/Alignments/RNASeq/"sorted_$(basename "$filename" .sam).bam"
done

for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_*
do
	samtools index $filename
done


for filename in $pathToParentDirectory/Termi/Alignments/TermSeq/trimmed_*.sam
do
	echo $filename 
	awk '/^@/ || $5 > "0"' $filename > $pathToParentDirectory/Termi/Alignments/TermSeq/"filtered_$(basename "$filename")"
done


for filename in $pathToParentDirectory/Termi/Alignments/TermSeq/filtered_*
do
	samtools view -bS $filename |
	samtools sort - -o $pathToParentDirectory/Termi/Alignments/TermSeq/"sorted_$(basename "$filename" .sam).bam"  > $pathToParentDirectory/Termi/Alignments/TermSeq/"sorted_$(basename "$filename" .sam).bam"
done

for filename in $pathToParentDirectory/Termi/Alignments/TermSeq/sorted_*
do
	samtools index $filename
done

printf "\n##########################################################"
printf '\n    all Files filtered and sorted'
printf "\n##########################################################\n\n"

# this is a good stopping point / checkpoint if you want to avoid re-doing huge computations while debugging.  The above code is the really computationally expensive stuff.  The remainder takes 1-2 hours.

###############################################################################################################################################

# calculate genome Coverage of RNA-Seq files and end-nuc counts of Term-Seq files

for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_*.bam
do
	genomeCoverageBed -d -split -ibam $filename > $pathToParentDirectory/Termi/Alignments/RNASeq/"$(basename $filename .bam)_genomeCoverageBed.bedgraph"
done

# add header with unique reads to bedgraph files

listOfNames=('ERX1320302_Bacillus_subtilis' 'ERR1248404_Enterococcus_faecalis_chromosome' 'ERR1248404_Enterococcus_faecalis_plasmid1' \
	'ERR1248404_Enterococcus_faecalis_plasmid2' 'ERR1248404_Enterococcus_faecalis_plasmid3' \
	'ERR1248439_Listeria_monocytogenes' \
	'316-RNASeqND0min-d-CACGAAT_Streptococcus_pneumoniae' '316-RNASeqND0min-c-ACAAGTT_Streptococcus_pneumoniae' \
	'316-RNASeqND0min-b-AATTCAT_Streptococcus_pneumoniae' '316-RNASeqND0min-a-AAGCAAT_Streptococcus_pneumoniae' )

for name in "${listOfNames[@]}"
do
	uniqueReads="$(awk '$5 > "0"' $pathToParentDirectory/Termi/Alignments/RNASeq/filtered_trimmed_${name}.sam | wc -l)"
	sed -i '1itrack type=bedGraph name=\"'$uniqueReads'\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20' \
	$pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${name}_genomeCoverageBed.bedgraph
done


$PYTHON $pathToParentDirectory/$SCRIPTS/5primePosfromSam_2_Bedgraph.py \
		-termSeq $pathToParentDirectory/Termi/Alignments/TermSeq/filtered_ \
		-o $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_


printf "\n##########################################################"
printf "\n   RNA-Seq coverage and Term-Seq 5\' end nucs counted"
printf "\n##########################################################\n\n"

################################################################################################################################################

mkdir -p $pathToParentDirectory"/Termi/Results"

# draw scatterplots of max Term-Seq end nuc counts vs. avg. RNA-Seq coverage over specified intervals of a set number of nucleotides 
# outputs Term-Seq counts overlapping genes (known class negative), counts overlapping known terminators (known class positive) and all counts

termSeqFiles1=('ERX1320300_Bacillus_subtilis' 'ERR1248401_Enterococcus_faecalis_chromosome' 'ERR1248401_Enterococcus_faecalis_plasmid1' \
	'ERR1248401_Enterococcus_faecalis_plasmid2' 'ERR1248401_Enterococcus_faecalis_plasmid3' \
	'ERR1248436_Listeria_monocytogenes')
termSeqFiles2=('ERX1320301_Bacillus_subtilis' 'ERR1248402_Enterococcus_faecalis_chromosome' 'ERR1248402_Enterococcus_faecalis_plasmid1' \
	'ERR1248402_Enterococcus_faecalis_plasmid2' 'ERR1248402_Enterococcus_faecalis_plasmid3' \
	'ERR1248437_Listeria_monocytogenes')
termSeqFiles3=('ERX1304415_Bacillus_subtilis' 'ERR1248403_Enterococcus_faecalis_chromosome' 'ERR1248403_Enterococcus_faecalis_plasmid1' \
	'ERR1248403_Enterococcus_faecalis_plasmid2' 'ERR1248403_Enterococcus_faecalis_plasmid3' \
	'ERR1248438_Listeria_monocytogenes')
RNAseqFiles=('ERX1320302_Bacillus_subtilis' 'ERR1248404_Enterococcus_faecalis_chromosome' 'ERR1248404_Enterococcus_faecalis_plasmid1' \
	'ERR1248404_Enterococcus_faecalis_plasmid2' 'ERR1248404_Enterococcus_faecalis_plasmid3' \
	'ERR1248439_Listeria_monocytogenes')
gffFiles=('Bacillus_subtilis.gff' 'Enterococcus_faecalis_chromosome.gff3' 'Enterococcus_faecalis_plasmid1.gff3'\
	'Enterococcus_faecalis_plasmid2.gff3' 'Enterococcus_faecalis_plasmid3.gff3'\
	'Listeria_monocytogenes.gff')
brev=('BS' 'EF_chrom' 'EF_pl1' 'EF_pl2' 'EF_pl3' 'LM')


########################

$PYTHON $pathToParentDirectory/$SCRIPTS/combine_RNASeq_SP.py \
		-i $pathToParentDirectory/Termi/Alignments/RNASeq/ \
		-o $pathToParentDirectory/Termi/Alignments/RNASeq/

$PYTHON $pathToParentDirectory/$SCRIPTS/ratioTermRNASeq_SP.py  \
		-gff $pathToParentDirectory/Termi/Genomes/Streptococcus_pneumoniae.gff3 \
		-ts $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_SRR7160964_Streptocoocus_pneumoniae.bedgraph \
			$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_SRR7160965_Streptocoocus_pneumoniae.bedgraph\
			$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_SRR7160966_Streptocoocus_pneumoniae.bedgraph \
			$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_SRR7160967_Streptocoocus_pneumoniae.bedgraph \
		-rs $pathToParentDirectory/Termi/Alignments/RNASeq/combined_316_RNASeqND0min_Streptococcus_pneumoniae_genomeCoverage.bedgraph \
		-bed $pathToParentDirectory/Termi/Alignments/knownTerminators_Bacillus_subtilis.bed \
		-c 100 \
		-o $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_SP

######################

# whole genome

for ((i=0;i<${#termSeqFiles1[@]};++i))
do
	$PYTHON $pathToParentDirectory/$SCRIPTS/ratioTermRNASeq.py  \
			-gff $pathToParentDirectory/Termi/Genomes/${gffFiles[i]} \
			-ts $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles1[i]}.bedgraph \
				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles2[i]}.bedgraph \
				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles3[i]}.bedgraph \
			-rs $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${RNAseqFiles[i]}_genomeCoverageBed.bedgraph \
			-bed $pathToParentDirectory/Termi/Alignments/knownTerminators_Bacillus_subtilis.bed \
			-c 100 \
			-o $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}
done

# first 500k and second 500k nucleotides for B.subtilis

# brev=('BS' 'EF_chrom' 'EF_pl1' 'EF_pl2' 'EF_pl3' 'LM')
# termSeqFiles1=('ERX1320300_Bacillus_subtilis')
# termSeqFiles2=('ERX1320301_Bacillus_subtilis')
# termSeqFiles3=('ERX1304415_Bacillus_subtilis')
# RNAseqFiles=('ERX1320302_Bacillus_subtilis')
# gffFiles=('Bacillus_subtilis.gff')
# brev=('BS')

# for ((i=0;i<${#termSeqFiles1[@]};++i))
# do
# 	$PYTHON $pathToParentDirectory/$SCRIPTS/ratioTermRNASeq.py  \
# 			-gff $pathToParentDirectory/Termi/Genomes/${gffFiles[i]} \
# 			-ts $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles1[i]}.bedgraph \
# 				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles2[i]}.bedgraph \
# 				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles3[i]}.bedgraph \
# 			-rs $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${RNAseqFiles[i]}_genomeCoverageBed.bedgraph \
# 			-bed $pathToParentDirectory/Termi/Alignments/knownTerminators_Bacillus_subtilis.bed \
# 			-c 100 \
# 			-g 1 500000\
# 			-o $pathToParentDirectory/Termi/Results/first500k_filtered_trim_scaled_${brev[i]}
# 	$PYTHON $pathToParentDirectory/$SCRIPTS/ratioTermRNASeq.py  \
# 			-gff $pathToParentDirectory/Termi/Genomes/${gffFiles[i]} \
# 			-ts $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles1[i]}.bedgraph \
# 				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles2[i]}.bedgraph \
# 				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles3[i]}.bedgraph \
# 			-rs $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${RNAseqFiles[i]}_genomeCoverageBed.bedgraph \
# 			-bed $pathToParentDirectory/Termi/Alignments/knownTerminators_Bacillus_subtilis.bed \
# 			-c 100 \
# 			-g 500000 1000000\
# 			-o $pathToParentDirectory/Termi/Results/second500k_filtered_trim_scaled_${brev[i]}
# done

printf "\n##########################################################"
printf '\n   max. Term-Seq vs. avg. RNA-Seq over all replicates'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# get fasta files from ratio results (all counts; 100 nucs long) from all species that are not B.subtilis
# run RNIE, filter out counts with RNIE scores over 20 as "known positives"

fastaFiles=('Enterococcus_faecalis_chromosome.fasta' 'Enterococcus_faecalis_plasmid1.fasta'\
	'Enterococcus_faecalis_plasmid2.fasta' 'Enterococcus_faecalis_plasmid3.fasta'\
	'Listeria_monocytogenes.fasta' 'Streptococcus_pneumoniae.fasta')
gffFiles2=('Enterococcus_faecalis_chromosome.gff3' 'Enterococcus_faecalis_plasmid1.gff3'\
	'Enterococcus_faecalis_plasmid2.gff3' 'Enterococcus_faecalis_plasmid3.gff3'\
	'Listeria_monocytogenes.gff' 'Streptococcus_pneumoniae.gff3')
brev2=('EF_chrom' 'EF_pl1' 'EF_pl2' 'EF_pl3' 'LM' 'SP')

for ((i=0;i<${#fastaFiles[@]};++i))
do
	echo $pathToParentDirectory/Termi/Genomes/${fastaFiles[i]}

	bedtools getfasta -fi $pathToParentDirectory/Termi/Genomes/${fastaFiles[i]} \
		-bed $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS.bed\
 		-name | \
	$PYTHON $pathToParentDirectory/$SCRIPTS/processBedtoolsGetFasta.py -revRev > $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS.fasta 
	$pathToParentDirectory/Termi/RNIE/rnie.pl \
		--gene \
		-f $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS.fasta \
	 	-md $pathToParentDirectory/Termi/RNIE/models \
		-th 0 \
	 	-p $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS_th0
	$PYTHON $pathToParentDirectory/$SCRIPTS/filterRNIE.py\
		-i $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff \
		-o $pathToParentDirectory/Termi/Results/wholeGenome_${brev2[i]}_
done


printf "\n##########################################################"
printf '\n            all RNIE scores filtered'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# classify counts into "predicted positives/terminators" and "predicted negatives"

mkdir -p $pathToParentDirectory"/Termi/Results/Classification"

# b.subtilis
$PYTHON $pathToParentDirectory/$SCRIPTS/classification.py  \
		-pos $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_BS_50_nucs_TScountsOverlappingTerminators \
		-neg $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_BS_50_nucs_TScountsOverlappingGenes \
		-all $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_BS_50_nucs_TSvsRS \
		-gff $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.gff  \
		-term $pathToParentDirectory/Termi/Alignments/knownTerminators_Bacillus_subtilis.bed \
		-o $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_

# rest
for ((i=0;i<${#gffFiles2[@]};++i))
do
	echo ${gffFiles2[i]} 
	$PYTHON $pathToParentDirectory/$SCRIPTS/classification.py  \
		-pos $pathToParentDirectory/Termi/Results/wholeGenome_${brev2[i]}_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed \
		-neg $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucs_TScountsOverlappingGenes \
		-all $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucs_TSvsRS \
		-gff $pathToParentDirectory/Termi/Genomes/${gffFiles2[i]} \
		-o $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_
done

printf "\n##########################################################"
printf '\n    classified predicted positives and negatives'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# # removing predictions over 150 nucleotides downstream of genes
brev=('BS' 'EF_chrom' 'EF_pl1' 'EF_pl2' 'EF_pl3' 'LM' 'SP')
gffFiles=('Bacillus_subtilis.gff' 'Enterococcus_faecalis_chromosome.gff3' 'Enterococcus_faecalis_plasmid1.gff3'\
	'Enterococcus_faecalis_plasmid2.gff3' 'Enterococcus_faecalis_plasmid3.gff3'\
	'Listeria_monocytogenes.gff' 'Streptococcus_pneumoniae.gff3')
mkdir -p $pathToParentDirectory"/Termi/Results/Distance"

for ((i=0;i<${#gffFiles[@]};++i))
do
	$PYTHON $pathToParentDirectory/$SCRIPTS/position.py \
		-pos $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_${brev[i]}_predictedTerminators_NO_knownTerminators_NO_genes_long.bed \
		-neg $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_${brev[i]}_predictedNegatives_NO_knownTerminators_NO_genes_long.bed \
		-gff $pathToParentDirectory/Termi/Genomes/${gffFiles[i]}  \
		-o $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_
done


printf "\n##########################################################"
printf '\n       distances to closest genes calculated'
printf "\n##########################################################\n\n"

##########################################################################################################################################

# removing predicted positives that are too similar to known terminators (BLAST bitscores > 30 )

# make blast db of known terminators
mkdir -p $pathToParentDirectory"/Termi/BLASTDB"

makeblastdb -in $pathToParentDirectory/Termi/Genomes/true_upper_UtoT_shortIDs.fa \
			-dbtype nucl -parse_seqids -out $pathToParentDirectory/Termi/BLASTDB/known_terminators_db

mkdir -p $pathToParentDirectory"/Termi/Results/BLAST"
mkdir -p $pathToParentDirectory"/Termi/Results/BLAST/withoutPolluted"


fastaFiles2=('Bacillus_subtilis.fasta' 'Enterococcus_faecalis_chromosome.fasta' 'Enterococcus_faecalis_plasmid1.fasta'\
	'Enterococcus_faecalis_plasmid2.fasta' 'Enterococcus_faecalis_plasmid3.fasta'\
	'Listeria_monocytogenes.fasta' 'Streptococcus_pneumoniae.fasta')

for ((i=0;i<${#fastaFiles2[@]};++i))
do
	bedtools getfasta \
		-fi $pathToParentDirectory/Termi/Genomes/${fastaFiles2[i]} \
		-bed $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_predictedTerminators_NO_knownTerminators_NO_genes_long.bed \
 		-name | \
	$PYTHON $pathToParentDirectory/$SCRIPTS/processBedtoolsGetFasta.py -revRev > $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta 
	blastn -query $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta \
		-db $pathToParentDirectory/Termi/BLASTDB/known_terminators_db \
		-word_size 7 \
		-outfmt 6 |
	sort -k12 -rn > $pathToParentDirectory/Termi/Results/BLAST/wholeGenome_filtered_trim_scaled_${brev[i]}_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.tab
	$PYTHON $pathToParentDirectory/$SCRIPTS/filterBLAST.py \
		-term  $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_predictedTerminators_NO_knownTerminators_NO_genes_long.bed \
		-blast $pathToParentDirectory/Termi/Results/BLAST/wholeGenome_filtered_trim_scaled_${brev[i]}_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.tab \
		-o $pathToParentDirectory/Termi/Results/BLAST/ 
done


for file in $pathToParentDirectory/Termi/Results/BLAST/BS*.bed
do
	sort -k2 -n $file > $pathToParentDirectory/Termi/Results/BLAST/sorted_"$(basename "$file")"
	awk '$2>1000000' $pathToParentDirectory/Termi/Results/BLAST/sorted_"$(basename "$file")" > $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_"$(basename "$file")"
done

for file in $pathToParentDirectory/Termi/Results/BLAST/{LM*,EF_chrom*}.bed
do
	sort -k2 -n $file > $pathToParentDirectory/Termi/Results/BLAST/sorted_"$(basename "$file")"
	awk '$2>500000' $pathToParentDirectory/Termi/Results/BLAST/sorted_"$(basename "$file")" > $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_"$(basename "$file")"
done

for file in $pathToParentDirectory/Termi/Results/BLAST/SP*.bed
do
	sort -k2 -n $file > $pathToParentDirectory/Termi/Results/BLAST/sorted_"$(basename "$file")"
	awk '$2>500000' $pathToParentDirectory/Termi/Results/BLAST/sorted_"$(basename "$file")" > $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_"$(basename "$file")"
done

printf "\n##########################################################"
printf '\n\t BLAST predicted against known Terminators'
printf "\n##########################################################\n\n"

##########################################################################################################################################

# using esl-shuffle with a 1st order Markov process to generate negatives and diresidue shuffle using Altschul-Erickson algorithm to generate sequences to embed negatives into

mkdir -p $pathToParentDirectory"/Termi/Results/Negatives"
mkdir -p $pathToParentDirectory"/Termi/Results/Negatives/1NegativePerPositive"
mkdir -p $pathToParentDirectory"/Termi/Results/Negatives/100NegativesPerPositive"


bedFiles=()


for file in $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/*.bed
do
	bedFiles=("${bedFiles[@]}" "$file")
done


fastaFiles3=('Bacillus_subtilis.fasta' 'Bacillus_subtilis.fasta' 'Bacillus_subtilis.fasta' \
	'Enterococcus_faecalis_chromosome.fasta' 'Enterococcus_faecalis_chromosome.fasta' 'Enterococcus_faecalis_chromosome.fasta' \
	'Listeria_monocytogenes.fasta' 'Listeria_monocytogenes.fasta' 'Listeria_monocytogenes.fasta' \
	'Streptococcus_pneumoniae.fasta' 'Streptococcus_pneumoniae.fasta' 'Streptococcus_pneumoniae.fasta')


 brev3=('BS' 'BS' 'BS' 'EF_chrom' 'EF_chrom' 'EF_chrom' 'LM' 'LM' 'LM' 'SP' 'SP' 'SP')

for ((i=0;i<${#bedFiles[@]};++i))
do
	bedtools getfasta \
	-fi $pathToParentDirectory/Termi/Genomes/${fastaFiles3[i]} \
	-bed ${bedFiles[i]} \
	-name | \
	$PYTHON $pathToParentDirectory/$SCRIPTS/processBedtoolsGetFasta.py -revRev > $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/"$(basename "${bedFiles[i]}" .bed)".fasta
done

# using esl-shuffle with a 1st order Markov process to generate negatives
for file in $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/*predTerm*.fasta
do
	esl-shuffle -N 1 -L 100 -1 --seed 1 $file > $pathToParentDirectory/Termi/Results/Negatives/1NegativePerPositive/1shuffled_"$(basename "$file")"
	esl-shuffle -N 100 -L 100 -1 --seed 1 $file > $pathToParentDirectory/Termi/Results/Negatives/100NegativesPerPositive/100shuffled_"$(basename "$file")"
done


# diresidue shuffle using Altschul-Erickson algorithm to generate sequences to embed negatives/predicted terminators into
for file in $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/*500*.fasta
do
	esl-shuffle -N 1 -L 500 -d --seed 1 $file > $pathToParentDirectory/Termi/Results/Negatives/1NegativePerPositive/1shuffled_"$(basename "$file")"
	esl-shuffle -N 100 -L 500 -d --seed 1 $file > $pathToParentDirectory/Termi/Results/Negatives/100NegativesPerPositive/100shuffled_"$(basename "$file")"
done

# removing linebreak in fasta files
for file in $pathToParentDirectory/Termi/Results/Negatives/1NegativePerPositive/*.fasta
do
	awk '!/^>/ { printf "%s", $0; n = "\n" } 
	/^>/ { print n $0; n = "" }
	END { printf "%s", n } 
	' $file \
	> $pathToParentDirectory/Termi/Results/Negatives/1NegativePerPositive/oneLine_"$(basename "$file")"
done

for file in $pathToParentDirectory/Termi/Results/Negatives/100NegativesPerPositive/*.fasta
do
	awk '!/^>/ { printf "%s", $0; n = "\n" } 
	/^>/ { print n $0; n = "" }
	END { printf "%s", n } 
	' $file \
	> $pathToParentDirectory/Termi/Results/Negatives/100NegativesPerPositive/oneLine_"$(basename "$file")"
done

mkdir -p $pathToParentDirectory"/Termi/Results/Embedded"
mkdir -p $pathToParentDirectory"/Termi/Results/Embedded/1NegativePerPositive"
mkdir -p $pathToParentDirectory"/Termi/Results/Embedded/100NegativesPerPositive"

# embed predicted terminators and shuffled negatives into permuted sequences
for ((i=0;i<${#brev3[@]};++i))
do
	$PYTHON $pathToParentDirectory/$SCRIPTS/embedding.py \
		-term $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_${brev3[i]}_predTerm_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta  \
		-neg $pathToParentDirectory/Termi/Results/Negatives/1NegativePerPositive/oneLine_1shuffled_cut_${brev3[i]}_predTerm_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta \
		-front $pathToParentDirectory/Termi/Results/Negatives/1NegativePerPositive/oneLine_1shuffled_cut_${brev3[i]}_500front_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta \
		-back $pathToParentDirectory/Termi/Results/Negatives/1NegativePerPositive/oneLine_1shuffled_cut_${brev3[i]}_500back_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta \
		-frontN $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_${brev3[i]}_500front_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta\
		-backN $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_${brev3[i]}_500back_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta\
		-o $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev3[i]}_

		$PYTHON $pathToParentDirectory/$SCRIPTS/embedding.py \
		-term $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_${brev3[i]}_predTerm_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta  \
		-neg $pathToParentDirectory/Termi/Results/Negatives/100NegativesPerPositive/oneLine_100shuffled_cut_${brev3[i]}_predTerm_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta \
		-front $pathToParentDirectory/Termi/Results/Negatives/100NegativesPerPositive/oneLine_100shuffled_cut_${brev3[i]}_500front_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta \
		-back $pathToParentDirectory/Termi/Results/Negatives/100NegativesPerPositive/oneLine_100shuffled_cut_${brev3[i]}_500back_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta \
		-frontN $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_${brev3[i]}_500front_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta\
		-backN $pathToParentDirectory/Termi/Results/BLAST/withoutPolluted/cut_${brev3[i]}_500back_BLAST_predictedTerminators_NO_knownTerminators_NO_genes_long.fasta\
		-o $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev3[i]}_
done



for file in $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/*.fasta
do
	cp $file  $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/"$(basename "$file" .fasta).csv"
done

for file in $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/*.fasta
do
	cp $file  $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/"$(basename "$file" .fasta).csv"
done

printf "\n##########################################################"
printf '\n\t sequences shuffled and embedded'
printf "\n##########################################################\n\n"

#############################################################################################################################################
