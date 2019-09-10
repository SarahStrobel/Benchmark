# !/bin/bash

# create new directory for project

pathToParentDirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

################################################################################################################################################

# download scripts from Sarah Strobel's github to new folder with wget 
wget https://raw.githubusercontent.com/SarahStrobel/Benchmark/master/Scripts/5primePosfromSam_2_Bedgraph.py -P $pathToParentDirectory"/Scripts/"
wget https://raw.githubusercontent.com/SarahStrobel/Benchmark/master/Scripts/classification.py -P $pathToParentDirectory"/Scripts/"
wget https://raw.githubusercontent.com/SarahStrobel/Benchmark/master/Scripts/position.py -P $pathToParentDirectory"/Scripts/"
wget https://raw.githubusercontent.com/SarahStrobel/Benchmark/master/Scripts/ratioTermRNASeq.py -P $pathToParentDirectory"/Scripts/"

printf "\n##########################################################"
printf '\n\t\t all scripts downloaded'
printf "\n##########################################################\n\n"

################################################################################################################################################

# download known terminator file from Paul Gardners's github to new folder with wget
wget https://raw.githubusercontent.com/ppgardne/RNIE-benchmark/master/training-test/true.fa -P $pathToParentDirectory"/knownTerminators/"

# change u to ts and make uppercase in known terminators sequences
sed 's/U/t/g' $pathToParentDirectory"/knownTerminators/true.fa" > $pathToParentDirectory"/knownTerminators/true_UtoT.fa" 
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1 " " $2}}' $pathToParentDirectory"/knownTerminators/true_UtoT.fa" \
	> $pathToParentDirectory"/knownTerminators/true_upper_UtoT.fa"

# download B.subtilis known terminator file form the Lin group homepage
wget http://lin-group.cn/server/iTerm-PseKNC/dependent_data2.csv -O $pathToParentDirectory"/knownTerminators/dependent_data2.fa"

# concat the terminator files
cat $pathToParentDirectory"/knownTerminators/dependent_data2.fa" $pathToParentDirectory"/knownTerminators/true_upper_UtoT.fa" \
	> $pathToParentDirectory"/knownTerminators/knownTerminators.fa"

printf "\n##########################################################"
printf '\n\t\t Known terminators downloaded'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# create new folders for RNA-Seq and Term-Seq data
mkdir -p $pathToParentDirectory"/RNASeq"
mkdir -p $pathToParentDirectory"/TermSeq"


# download RNA-Seq and Term-Seq data with sra toolkits fastq-dump
fastq-dump ERX1320302 --outdir $pathToParentDirectory"/RNASeq"
fastq-dump ERX1304415 --outdir $pathToParentDirectory"/TermSeq"
fastq-dump ERX1320300 --outdir $pathToParentDirectory"/TermSeq"
fastq-dump ERX1320301 --outdir $pathToParentDirectory"/TermSeq"

printf "\n##########################################################"
printf '\n    all RNA-Seq and Term-Seq experiments downloaded'
printf "\n##########################################################\n\n"

# ###############################################################################################################################################

# download B.subtilis genome and gene annotation file from Sarah Strobel's github to new folder with wget
wget https://raw.githubusercontent.com/SarahStrobel/Genomes/master/Bacillus_subtilis_UtoT.fasta -P $pathToParentDirectory"/Genomes"
wget https://raw.githubusercontent.com/SarahStrobel/Genomes/master/GCF_000009045.1_ASM904v1_genomic.gff -P $pathToParentDirectory"/Genomes"

printf "\n##########################################################"
printf '\n\t\t B.subtilis genome downloaded'
printf "\n##########################################################\n\n"

###############################################################################################################################################

mkdir -p $pathToParentDirectory"/Alignments"
mkdir -p $pathToParentDirectory"/Alignments/TermSeq"
mkdir -p $pathToParentDirectory"/Alignments/RNASeq"

# creating Index files of B.subtilis genome with Novoindex for alignments with Novoalign
novoindex $pathToParentDirectory"/Alignments/Bacillus_subtilis_UtoT_index" $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta"

# align known Terminators to B.subtilis genomes with Novoalign
novoalign -f $pathToParentDirectory"/knownTerminators/knownTerminators.fa" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/Bacillus_subtilis_knownTerminators.sam"
printf "\n"

# align Term-Seq and RNA-Seq files to B.subtilis genome with Novoalign
novoalign -f $pathToParentDirectory"/RNASeq/ERX1320302.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.sam"
printf "\n"
novoalign -f $pathToParentDirectory"/TermSeq/ERX1304415.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.sam"
printf "\n"
novoalign -f $pathToParentDirectory"/TermSeq/ERX1320300.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/TermSeq/ERX1320300.sam"
printf "\n"
novoalign -f $pathToParentDirectory"/TermSeq/ERX1320301.fastq" -d $pathToParentDirectory/"Alignments/Bacillus_subtilis_UtoT_index" \
			-o SAM > $pathToParentDirectory"/Alignments/TermSeq/ERX1320301.sam"

printf "\n##########################################################"
printf '\n\t all files aligned to B.subitils genome'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# converting into machine readable bam format with samtools view
# # one line gives error --> delete
sed -i '428d' $pathToParentDirectory"/Alignments/Bacillus_subtilis_knownTerminators.sam"
samtools view -S -b $pathToParentDirectory"/Alignments/Bacillus_subtilis_knownTerminators.sam" > $pathToParentDirectory/"Alignments/Bacillus_subtilis_knownTerminators.bam"

samtools view -S -b $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.sam" > $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.bam"
samtools view -S -b $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.sam" > $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.bam"
samtools view -S -b $pathToParentDirectory"/Alignments/TermSeq/ERX1320300.sam" > $pathToParentDirectory"/Alignments/TermSeq/ERX1320300.bam"
samtools view -S -b $pathToParentDirectory"/Alignments/TermSeq/ERX1320301.sam" > $pathToParentDirectory"/Alignments/TermSeq/ERX1320301.bam"

# sorting and indexing bam files to make them readable for genome browsers (e.g. IGV)
samtools sort $pathToParentDirectory"/Alignments/RNASeq/ERX1320302.bam" -o $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam" \
				> $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam"
samtools sort $pathToParentDirectory"/Alignments/TermSeq/ERX1304415.bam" -o $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1304415.bam" \
				> $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1304415.bam"
samtools sort $pathToParentDirectory"/Alignments/TermSeq/ERX1320300.bam" -o $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320300.bam" \
				> $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320300.bam"
samtools sort $pathToParentDirectory"/Alignments/TermSeq/ERX1320301.bam" -o $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320301.bam" \
				> $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320301.bam"

samtools index $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam"
samtools index $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1304415.bam"
samtools index $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320300.bam"
samtools index $pathToParentDirectory"/Alignments/TermSeq/sorted_ERX1320301.bam"

printf "\n##########################################################"
printf '\n\t\t all files sorted and indexed'
printf "\n##########################################################\n"

###############################################################################################################################################

# convert known Terminator bam file to bed file 
bedtools bamtobed -i $pathToParentDirectory"/Alignments/Bacillus_subtilis_knownTerminators.bam" \
					> $pathToParentDirectory"/Alignments/Bacillus_subtilis_knownTerminators.bed"

###############################################################################################################################################

# calculate genome coverage of RNA-Seq file with 
genomeCoverageBed -d -split -ibam $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302.bam" \
					> $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302_genomeCoverageBed.bedgraph"

##############################################################################################################################################

# calculate 5' end nuc count from Term-Seq alignment files, can be viewed in geneome browser (e.g. IGV)
python $pathToParentDirectory"/Scripts/5primePosfromSam_2_Bedgraph.py" -rnaSeq $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302_genomeCoverageBed.bedgraph" \
										-termSeq $pathToParentDirectory"/Alignments/TermSeq/" \
										-o $pathToParentDirectory"/Alignments/TermSeq/stop_"

printf "\n##########################################################"
printf "\n   RNA-Seq coverage and Term-Seq 5\' end nucs counted"
printf "\n##########################################################\n\n"

###############################################################################################################################################

mkdir $pathToParentDirectory"/Results"

# draw scatterplots of max Term-Seq end nuc counts vs. avg. RNA-Seq coverage over specified intervals of a set number of nucleotides 
# outputs Term-Seq counts overlapping genes (known class negative) and overlapping known terminators (known class positive) and all points

# required:
		# - Term-Seq end nuc count replicates (bedgraph files)
		# - RNA-Seq coverage file (bedgraph files)
		# - known Terminator (bed file)
		# - gene annotation (gff file)
		
# optional: 
		# - split Genome into chunks (e.g. only first 500k nucs)
		# - number of nucleotides to chop off both sides of genes (default 100)
		# - decide the lenght of intervals to cut the genome into (default 50)
		# - decide if you want to calculate the max or average in the RNA-Seq data (default avgerage)
		# - decide over how many nucs you want to average/max the RNA-Seq data over (default 50)

# example for whole genome in B.subtilis, 50 nuc intervals, rna seq average in span of 50 nucs 
python $pathToParentDirectory"/Scripts/ratioTermRNASeq.py" -gff $pathToParentDirectory"/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
						-ts $pathToParentDirectory"/Alignments/TermSeq/stop_ERX1320300.bedgraph" \
							$pathToParentDirectory"/Alignments/TermSeq/stop_ERX1320301.bedgraph" \
							$pathToParentDirectory"/Alignments/TermSeq/stop_ERX1304415.bedgraph" \
						-rs $pathToParentDirectory"/Alignments/RNASeq/sorted_ERX1320302_genomeCoverageBed.bedgraph" \
						-bed $pathToParentDirectory"/Alignments/Bacillus_subtilis_knownTerminators.bed" \
						-g 1 4215606 \
						-o $pathToParentDirectory"/Results/wholeGenome"

printf "\n##########################################################"
printf '\n   max. Term-Seq vs. avg. RNA-Seq over all replicates'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# drawing decision boundary for classification into "predicted terminators" and "predicted negatives"
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

printf "\n##########################################################"
printf '\n\t\t decision boundary drawn'
printf "\n##########################################################\n\n"

###############################################################################################################################################

# removing predictions over 150 nucleotides downstream of genes 
# outputs scatterplots with color coded distances 
mkdir $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/"

python $pathToParentDirectory"/Scripts/position.py" \
		-pos $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
		-neg $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
		-gene $pathToParentDirectory"/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
		-o $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/"

printf "\n##########################################################"
printf '\n\t distance to closest genes calculated'
printf "\n##########################################################\n\n"

#############################################################################################################################################

# bedtools get fasta to get sequences and making reverse complements for - strand
bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
 		-name | \
cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		

bedtools getfasta -fi $pathToParentDirectory"/Genomes/Bacillus_subtilis_UtoT.fasta" \
		-bed $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
 		-name | \
cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedNegatives_NO_knownTerminators_NO_genes.fasta" 


printf "\n##########################################################"
printf '\n\t\t fasta files generated'
printf "\n##########################################################\n\n"

################################################################################################################################################

# blast predicted Terminators against known Terminators and sort by bitscore
blastn -query $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
		-db ~/blastdb/new_known_terminators_db \
		-word_size 7 \
		-outfmt 6 |
sort -k12 -rn > $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/BLAST/BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"

printf "\n##########################################################"
printf '\n\t BLAST predicted against known Terminators'
printf "\n##########################################################\n\n"

################################################################################################################################################

# remove predictions with BLAST bitscores over 30 
python $pathToParentDirectory"/Scripts/filterBLAST.py" -term  $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/Distance/Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
		-blast $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/BLAST/BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
		-o $pathToParentDirectory"/Results/PredictedPositivesAndNegatives/BLAST/"