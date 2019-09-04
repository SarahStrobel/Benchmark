# !/bin/bash

# # create new directory for project

pathToParentDirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

################################################################################################################################################

wget https://github.com/SarahStrobel/Benchmark/blob/master/Scripts/5primePosfromSam_2_Bedgraph.py -P $pathToParentDirectory"/Scripts/"
wget https://github.com/SarahStrobel/Benchmark/blob/master/Scripts/classification.py -P $pathToParentDirectory"/Scripts/"
wget https://github.com/SarahStrobel/Benchmark/blob/master/Scripts/position.py -P $pathToParentDirectory"/Scripts/"
wget https://github.com/SarahStrobel/Benchmark/blob/master/Scripts/ratioTermRNASeq.py -P $pathToParentDirectory"/Scripts/"

################################################################################################################################################

# download known terminator file from Paul Gardners's github to new folder with wget
wget https://raw.githubusercontent.com/ppgardne/RNIE-benchmark/master/training-test/true.fa -P $pathToParentDirectory"/knownTerminators/"

# change u to ts and make uppercase in known terminators sequences
sed 's/U/t/g' $pathToParentDirectory$knownTerminators"/true.fa" > $pathToParentDirectory"/knownTerminators/true_UtoT.fa" 
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1 " " $2}}' $pathToParentDirectory"/knownTerminators/true_UtoT.fa" \
	> $pathToParentDirectory"/knownTerminators/true_upper_UtoT.fa"

###############################################################################################################################################

# create new folders for RNA-Seq and Term-Seq data
mkdir -p $pathToParentDirectory"/RNASeq"
mkdir -p $pathToParentDirectory"/TermSeq"


# download RNA-Seq and Term-Seq data with sra toolkits fastq-dump
fastq-dump ERX1320302 --outdir $pathToParentDirectory"/RNASeq"
fastq-dump ERX1304415 --outdir $pathToParentDirectory"/TermSeq"
fastq-dump ERX1320300 --outdir $pathToParentDirectory"/TermSeq"
fastq-dump ERX1320301 --outdir $pathToParentDirectory"/TermSeq"

###############################################################################################################################################

# download B.subtilis and E.coli genomes from Sarah Strobel's github to new folder with wget
wget https://raw.githubusercontent.com/SarahStrobel/Genomes/master/Bacillus_subtilis_UtoT.fasta -P $pathToParentDirectory"/Genomes"
wget https://raw.githubusercontent.com/SarahStrobel/Genomes/master/Escherichia_coli_IAI39.fasta -P $pathToParentDirectory"/Genomes"
wget https://raw.githubusercontent.com/SarahStrobel/Genomes/master/GCF_000009045.1_ASM904v1_genomic.gff -P $pathToParentDirectory"/Genomes"

###############################################################################################################################################

mkdir -p $pathToParentDirectory"/Alignments"
mkdir -p $pathToParentDirectory"/Alignments/TermSeq"
mkdir -p $pathToParentDirectory"/Alignments/RNASeq"


# creating Index files with Novoindex for Alignments with Novoalign
novoindex $pathToParentDirectory"/Alignments/Bacillus_subtilis_UtoT_index" $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta"
novoindex $pathToParentDirectory"/Alignments/Escherichia_coli_index" $pathToParentDirectory"/Genomes/Escherichia_coli_IAI39.fasta"


# align known Terminators to B.subtilis and E.coli genomes with Novoalign
novoalign -f $pathToParentDirectory"/knownTerminators/true_upper_UtoT.fa" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/Bacillus_subtilis_true_upper_UtoT.sam"
novoalign -f $pathToParentDirectory"/knownTerminators/true_upper_UtoT.fa" -d $pathToParentDirectory/"Alignments/Escherichia_coli_index" \
			-o SAM > $pathToParentDirectory"/Alignments/Escherichia_coli_true_upper_UtoT.sam"


# align Term-Seq and RNA-Seq files to B.subtilis genome with Novoalign
novoalign -f $pathToParentDirectory"/RNASeq/ERX1320302.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.sam"
novoalign -f $pathToParentDirectory"/TermSeq/ERX1304415.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.sam"
novoalign -f $pathToParentDirectory"/TermSeq/ERX1320300.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/TermSeq/ERX1320300.sam"
novoalign -f $pathToParentDirectory"/TermSeq/ERX1320301.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/TermSeq/ERX1320301.sam"


###############################################################################################################################################

# converting into machine readable bam format with samtools view
samtools view -S -b $pathToParentDirectory"/Alignments/Bacillus_subtilis_true_upper_UtoT.sam" \
				> $pathToParentDirectory"/Alignments/Bacillus_subtilis_true_upper_UtoT.bam"
samtools view -S -b $pathToParentDirectory"/Alignments/Escherichia_coli_true_upper_UtoT.sam" \
				> $pathToParentDirectory"/Alignments/Escherichia_coli_true_upper_UtoT.bam"

samtools view -S -b $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.sam" > $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.bam"
samtools view -S -b $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.sam" > $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.bam"
samtools view -S -b $pathToParentDirectory"/Alignments/TermSeq/ERX1320300.sam" > $pathToParentDirectory"/Alignments/TermSeq/ERX1320300.bam"
samtools view -S -b $pathToParentDirectory"/Alignments/TermSeq/ERX1320301.sam" > $pathToParentDirectory"/Alignments/TermSeq/ERX1320301.bam"



# sorting and indexing bam files to make them readable for igv
samtools sort $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.bam" -o $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam" \
				> $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam"
samtools sort $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.bam" -o $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1304415.bam" \
				> $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1304415.bam"
samtools sort $pathToParentDirectory"/Alignments/RNASeq/ERX1320300.bam" -o $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320300.bam" \
				> $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320300.bam"
samtools sort $pathToParentDirectory"/Alignments/RNASeq/ERX1320301.bam" -o $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320301.bam" \
				> $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320301.bam"


samtools index $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam"
samtools index $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1304415.bam"
samtools index $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320300.bam"
samtools index $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320301.bam"

###############################################################################################################################################

# convert known Terminator bam file to bed file 
bedtools bamtobed -i /$pathToParentDirectory"/Alignments/Bacillus_subtilis_true_upper_UtoT.bam" \
					> $pathToParentDirectory"/Alignments/Bacillus_subtilis_true_upper_UtoT.bed"

###############################################################################################################################################

# calculate genome coverage of RNA-Seq file with 
genomeCoverageBed -d -split -ibam $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam" \
					> $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302_genomeCoverageBed.bedgraph"

###############################################################################################################################################

# calculate 5' end nuc count from Term-Seq alignment files
python $pathToParentDirectory"/Scripts/5primePosfromSam_2_Bedgraph.py" -rnaSeq $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302_genomeCoverageBed.bedgraph" \
										-termSeq $pathToParentDirectory"/Alignments/TermSeq/" -o $pathToParentDirectory"/Alignments/TermSeq/stop_"

python $pathToParentDirectory"/Scripts/5primePosfromSam_2_Bedgraph.py" -rnaSeq "scr/k70san3/stsarah/Novoalign_Alignments/sorted_ERX1320302_genomeCoverage_d_UtoT.bedgraph" \
										-termSeq "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/" -o $pathToParentDirectory"/Alignments/TermSeq/stop_"	

###############################################################################################################################################

mkdir $pathToParentDirectory"/Results"

###############################################################################################################################################

# draw scatterplots of max Term-Seq end nuc counts vs. avg. RNA-Seq coverage over specified intervals of a set number of nucleotides 
# outputs Term-Seq counts overlapping genes (known class negative) and overlapping known terminators (known class positive) and all points

# takes:
		# - Term-Seq end nuc count replicates (bedgraph files)
		# - RNA-Seq coverage file (bedgraph files)
		# - known Terminator (bed file)
		# - gene annotation (gff file)
		

# optional: 
		# - split Genome into chunks (e.g. only first 500k nucs)
		# - decide the lenght of intervals to cut the genome into (default 50)
		# - decide if you want to calculate the max or average in the RNA-Seq data (default max)
		# - decide over how many nucs you want to average/max the RNA-Seq data over (default 50)

# example for first 500k nucleotides in B.subtilis
# python $pathToParentDirectory"/Scripts/ratioTermRNASeq.py" -gff $pathToParentDirectory"/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
# 						-ts $pathToParentDirectory"/Alignments/TermSeq/stop_ERX1320300.bedgraph" \
# 							$pathToParentDirectory"/Alignments/TermSeq/stop_ERX1320301.bedgraph" \
# 							$pathToParentDirectory"/Alignments/TermSeq/stop_ERX1304415.bedgraph" \
# 						-rs $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/Bacillus_subtilis_true_upper_UtoT.bed" \
# 						-g 1 500000 \
# 						-o  $pathToParentDirectory"/Results/first500k"

# example for whole genome in B.subtilis
python $pathToParentDirectory"/Scripts/ratioTermRNASeq.py" -gff $pathToParentDirectory"/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
						-ts $pathToParentDirectory"/Alignments/TermSeq/stop_ERX1320300.bedgraph" \
							$pathToParentDirectory"/Alignments/TermSeq/stop_ERX1320301.bedgraph" \
							$pathToParentDirectory"/Alignments/TermSeq/stop_ERX1304415.bedgraph" \
						-rs $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302_genomeCoverageBed.bedgraph" \
						-o $pathToParentDirectory"/Results/wholeGenome"


###############################################################################################################################################

# drawing decision boundary for classification into predicted terminators and predicted negatives
# outputs:
		#- PPV, TPR, FPR, TP, FP, TN, FN and overlaps on the console
		#- scatterplots with decision boundary
		#- bed files for predicted negatives and positives (not overlapping genes or not overlapping genes and not overlapping terminators)
		#- bed files creating 120 and 60 nucleotide long sequences (100/50 upstream and 20/10 downstream) for further investigation (e.g. BLAST, RNIE, mfold,...)

mkdir $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/"

python $pathToParentDirectory"/Scripts/classification.py" -pos $pathToParentDirectory"/Results/wholeGenome_50_nucs_TScountsOverlappingTerminators" \
						-neg $pathToParentDirectory"/Results/wholeGenome_50_nucs_TScountsOverlappingGenes" \
						-all $pathToParentDirectory"/Results/wholeGenome_50_nucs_TSvsRS" \
						-gene $pathToParentDirectory"/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
						-term $pathToParentDirectory"/Alignments/Bacillus_subtilis_true_upper_UtoT.bed" \
						-o $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/"


###############################################################################################################################################

# # bedtools get fasta to get sequences and making reverse complements for - strand

# bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
# 		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/120_predictedTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/120_predictedTerminators_NO_genes.fasta" 		

# bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
# 		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/120_predictedNegatives_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/120_predictedNegatives_NO_genes.fasta" 		

# bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
# 		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedTerminators_NO_genes.fasta" 		

# bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
# 		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedNegatives_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedNegatives_NO_genes.fasta" 	

###############################################################################################################################################

# looking for predictions up to 150 nucleotides downstream of genes 

mkdir $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/"

python $pathToParentDirectory"/Scripts/position.py" \
		-pos $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedTerminators_NO_genes.bed"\
		-neg $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedNegatives_NO_genes.bed"\
		-gene $pathToParentDirectory"/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
		-o $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/"

###############################################################################################################################################

# bedtools get fasta to get sequences and making reverse complements for - strand

bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_120_predictedTerminators_NO_genes.bed"\
 		-name | \
cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_120_predictedTerminators_NO_genes.fasta" 		

bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_120_predictedNegatives_NO_genes.bed"\
 		-name | \
cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_120_predictedNegatives_NO_genes.fasta"

bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedTerminators_NO_genes.bed"\
 		-name | \
cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedTerminators_NO_genes.fasta" 		

bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedNegatives_NO_genes.bed"\
 		-name | \
cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedNegatives_NO_genes.fasta" 


###############################################################################################################################################