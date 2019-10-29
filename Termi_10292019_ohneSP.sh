# !/bin/bash

pathToParentDirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
pathToRNIE="/homes/brauerei/stsarah/Masterarbeit/Infernal/"
pathToBLASTDB="/homes/brauerei/stsarah/blastdb/"

#################################################################################################################################################

# # download scripts, known terminators and genomes from Sarah Strobel's github to new folder 
# git clone https://github.com/SarahStrobel/Benchmark/ $pathToParentDirectory"/Termi/"

# printf "\n##########################################################"
# printf '\n all scripts, terminators and genomes downloaded'
# printf "\n##########################################################\n\n"

# #################################################################################################################################################

# # # download known terminator file from Paul Gardners's github to new folder with wget
# wget https://raw.githubusercontent.com/ppgardne/RNIE-benchmark/master/training-test/true.fa -P $pathToParentDirectory"/Termi/knownTerminators/"

# # change u to ts and make uppercase in known terminators sequences
# sed 's/U/t/g' $pathToParentDirectory"/Termi/knownTerminators/true.fa" |
# awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1 " " $2}}'  \
# 	> $pathToParentDirectory"/Termi/knownTerminators/true_upper_UtoT.fa"

# # download B.subtilis known terminator file form the Lin group homepage
# wget http://lin-group.cn/server/iTerm-PseKNC/dependent_data2.csv -O $pathToParentDirectory"/Termi/knownTerminators/dependent_data2.fa"

# printf "\n##########################################################"
# printf '\n\t\t Known terminators downloaded'
# printf "\n##########################################################\n\n"


# ################################################################################################################################################

# # # download fastq files wit fasterq-dump
# mkdir -p $pathToParentDirectory"/Termi/RNASeq"
# mkdir -p $pathToParentDirectory"/Termi/TermSeq"

# TermseqFiles=( ERX1304415 ERX1320300 ERX1320301 ERR1248401 ERR1248402 ERR1248403 ERR1248436 ERR1248437 ERR1248438 )
# for i in "${TermseqFiles[@]}"
# do
# 	echo $i
# 	fasterq-dump $i --outdir $pathToParentDirectory"/Termi/TermSeq"
# done

# RNAseqFiles=( ERX1320302 ERR1248404 ERR1248439 )

# for i in "${RNAseqFiles[@]}"
# do
# 	echo $i
# 	fasterq-dump $i --outdir $pathToParentDirectory"/Termi/RNASeq"
# done


# printf "\n##########################################################"
# printf '\n    all RNA-Seq and Term-Seq experiments downloaded'
# printf "\n##########################################################\n\n"

################################################################################################################################################

# # trim fastq files (quality control and adapter trimming)

# for filename in $pathToParentDirectory/Termi/RNASeq/*.fastq
# do
# 	fastp -i $filename -o  $pathToParentDirectory/Termi/RNASeq/"trimmed_$(basename "$filename" )"
# done

# for filename in $pathToParentDirectory/Termi/TermSeq/*.fastq
# do
# 	fastp -i $filename -o $pathToParentDirectory/Termi/TermSeq/"trimmed_$(basename "$filename" )"
# done

# printf "\n##########################################################"
# printf '\n    Adapter Sequences trimmed'
# printf "\n##########################################################\n\n"

################################################################################################################################################

# create index files and align fastq files to respective genomes with novoalign

# for filename in $pathToParentDirectory/Termi/Genomes/*.fasta
# do
# 	novoindex $pathToParentDirectory/Termi/Genomes/"$(basename "$filename" .fasta).index" $filename
# done

# mkdir $pathToParentDirectory/Termi/Alignments/
# mkdir $pathToParentDirectory/Termi/Alignments/RNASeq/
# mkdir $pathToParentDirectory/Termi/Alignments/TermSeq/

# # rna seq

# novoalign -f $pathToParentDirectory/Termi/RNASeq/trimmed_ERX1320302.fastq \
# 			-d $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.index \
# 			-o SAM > $pathToParentDirectory/Termi/Alignments/RNASeq/trimmed_ERX1320302_Bacillus_subtilis.sam

# novoalign -f $pathToParentDirectory/Termi/RNASeq/trimmed_ERR1248439.fastq \
# 			-d $pathToParentDirectory/Termi/Genomes/Listeria_monocytogenes.index \
# 			-o SAM > $pathToParentDirectory/Termi/Alignments/RNASeq/trimmed_ERR1248439_Listeria_monocytogenes.sam

# for filename in $pathToParentDirectory/Termi/Genomes/Enterococcus*.index
# do
# 	echo $pathToParentDirectory/Termi/RNASeq/trimmed_ERR1248404_"$(basename "$filename" .index).sam"
# 	novoalign 	-f $pathToParentDirectory/Termi/RNASeq/trimmed_ERR1248404.fastq \
# 				-d $filename \
# 				-o SAM > $pathToParentDirectory/Termi/Alignments/RNASeq/trimmed_ERR1248404_"$(basename "$filename" .index).sam"
# done



# # term seq

# for filename in $pathToParentDirectory/Termi/TermSeq/trimmed_ERX*
# do
# 	echo $pathToParentDirectory/Termi/RNASeq/"$(basename $filename .fastq)_Bacillus_subtilis.sam"
# 	novoalign -f $filename \
# 			-d $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.index \
# 			-o SAM > $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename $filename .fastq)_Bacillus_subtilis.sam"
# done

# for filename in $pathToParentDirectory/Termi/TermSeq/trimmed_ERR124843*
# do
# 	# echo $pathToParentDirectory/Termi/RNASeq/"$(basename $filename .fastq)_Listeria_monocytogenes.sam"
# 	novoalign -f $filename \
# 			-d $pathToParentDirectory/Termi/Genomes/Listeria_monocytogenes.index \
# 			-o SAM > $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename $filename .fastq)_Listeria_monocytogenes.sam"
# done

# for filename in $pathToParentDirectory/Termi/Genomes/Enterococcus*.index
# do
# 	for filename2 in $pathToParentDirectory/Termi/TermSeq/trimmed_ERR124840*
# 	do
# 		# echo $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename "$filename2" .fastq)"_$(basename "$filename" .index)".sam"
# 		novoalign 	-f $filename2 \
# 					-d $filename \
# 					-o SAM > $pathToParentDirectory/Termi/Alignments/TermSeq/"$(basename "$filename2" .fastq)"_$(basename "$filename" .index)".sam"
# 	done
# done



# # known terminators to B.subtilis

# novoalign 	-f $pathToParentDirectory/Termi/knownTerminators/true_upper_UtoT.fa \
# 			-d $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.index \
# 			-o SAM > $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.sam"


# samtools view -bS $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.sam" > $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.bam"

# bedtools bamtobed -i $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.bam" \
# 					> $pathToParentDirectory/Termi/Alignments/"knownTerminators_Bacillus_subtilis.bed"

# printf "\n##########################################################"
# printf '\n    all Files aligned'
# printf "\n##########################################################\n\n"

################################################################################################################################################

# # remove multimapped reads from sam files and use samtools view to convert into bam files, sort and index (viewable in igv)

# for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/*
# do
# 	echo $filename 
# 	awk '/^@/ || $5 > "0"' $filename > $pathToParentDirectory/Termi/Alignments/RNASeq/"filtered_$(basename "$filename")"
# done


# for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/filtered_*
# do
# 	samtools view -bS $filename |
# 	samtools sort - -o $pathToParentDirectory/Termi/Alignments/RNASeq/"sorted_$(basename "$filename" .sam).bam"  > $pathToParentDirectory/Termi/Alignments/RNASeq/"sorted_$(basename "$filename" .sam).bam"
# done

# for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_*
# do
# 	samtools index $filename
# done


# for filename in $pathToParentDirectory/Termi/Alignments/TermSeq/*
# do
# 	echo $filename 
# 	awk '/^@/ || $5 > "0"' $filename > $pathToParentDirectory/Termi/Alignments/TermSeq/"filtered_$(basename "$filename")"
# done


# for filename in $pathToParentDirectory/Termi/Alignments/TermSeq/filtered_*
# do
# 	samtools view -bS $filename |
# 	samtools sort - -o $pathToParentDirectory/Termi/Alignments/TermSeq/"sorted_$(basename "$filename" .sam).bam"  > $pathToParentDirectory/Termi/Alignments/TermSeq/"sorted_$(basename "$filename" .sam).bam"
# done

# for filename in $pathToParentDirectory/Termi/Alignments/TermSeq/sorted_*
# do
# 	samtools index $filename
# done

# printf "\n##########################################################"
# printf '\n    all Files filtered and sorted'
# printf "\n##########################################################\n\n"

################################################################################################################################################

# # calculate genome Coverage of RNA-Seq files and end-nuc counts of Term-Seq files

# for filename in $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_*.bam
# do
# 	echo $pathToParentDirectory/Termi/Alignments/RNASeq/"$(basename $filename .bam)_genomeCoverageBed.bedgraph"
# 	genomeCoverageBed -d -split -ibam $filename > $pathToParentDirectory/Termi/Alignments/RNASeq/"$(basename $filename .bam)_genomeCoverageBed.bedgraph"
# done


# listOfNames=('ERR1248404_Enterococcus_faecalis_chromosome' 'ERR1248404_Enterococcus_faecalis_plasmid1' \
# 	'ERR1248404_Enterococcus_faecalis_plasmid2' 'ERR1248404_Enterococcus_faecalis_plasmid3' \
# 	'ERR1248439_Listeria_monocytogenes' 'ERX1320302_Bacillus_subtilis')

# for name in "${listOfNames[@]}"
# do
# 	# echo $name
# 	# echo $pathToParentDirectory/Termi/Alignments/RNASeq/filtered_trimmed_${name}.sam
# 	uniqueReads="$(awk '$5 > "0"' $pathToParentDirectory/Termi/Alignments/RNASeq/filtered_trimmed_${name}.sam | wc -l)"
# 	echo $uniqueReads
# 	echo $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${name}_genomeCoverageBed.bedgraph
# 	sed -i '1itrack type=bedGraph name=\"'$uniqueReads'\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20' \
# 	$pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${name}_genomeCoverageBed.bedgraph
# done


# python $pathToParentDirectory/Termi/Scripts/5primePosfromSam_2_Bedgraph.py \
# 		-termSeq $pathToParentDirectory/Termi/Alignments/TermSeq/filtered_ \
# 		-o $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_


# printf "\n##########################################################"
# printf "\n   RNA-Seq coverage and Term-Seq 5\' end nucs counted"
# printf "\n##########################################################\n\n"

################################################################################################################################################

# mkdir $pathToParentDirectory"/Termi/Results"

# # draw scatterplots of max Term-Seq end nuc counts vs. avg. RNA-Seq coverage over specified intervals of a set number of nucleotides 
# # outputs Term-Seq counts overlapping genes (known class negative), counts overlapping known terminators (known class positive) and all counts

# termSeqFiles1=('ERR1248401_Enterococcus_faecalis_chromosome' 'ERR1248401_Enterococcus_faecalis_plasmid1' \
# 	'ERR1248401_Enterococcus_faecalis_plasmid2' 'ERR1248401_Enterococcus_faecalis_plasmid3' \
# 	'ERR1248436_Listeria_monocytogenes' 'ERX1320300_Bacillus_subtilis')
# termSeqFiles2=('ERR1248402_Enterococcus_faecalis_chromosome' 'ERR1248402_Enterococcus_faecalis_plasmid1' \
# 	'ERR1248402_Enterococcus_faecalis_plasmid2' 'ERR1248402_Enterococcus_faecalis_plasmid3' \
# 	'ERR1248437_Listeria_monocytogenes' 'ERX1320301_Bacillus_subtilis')
# termSeqFiles3=('ERR1248403_Enterococcus_faecalis_chromosome' 'ERR1248403_Enterococcus_faecalis_plasmid1' \
# 	'ERR1248403_Enterococcus_faecalis_plasmid2' 'ERR1248403_Enterococcus_faecalis_plasmid3' \
# 	'ERR1248438_Listeria_monocytogenes' 'ERX1304415_Bacillus_subtilis')
# RNAseqFiles=('ERR1248404_Enterococcus_faecalis_chromosome' 'ERR1248404_Enterococcus_faecalis_plasmid1' \
# 	'ERR1248404_Enterococcus_faecalis_plasmid2' 'ERR1248404_Enterococcus_faecalis_plasmid3' \
# 	'ERR1248439_Listeria_monocytogenes' 'ERX1320302_Bacillus_subtilis')
gffFiles=('Enterococcus_faecalis_Chromosome.gff3' 'Enterococcus_faecalis_plasmid1.gff3'\
	'Enterococcus_faecalis_plasmid2.gff3' 'Enterococcus_faecalis_plasmid3.gff3'\
	'Listeria_monocytogenes.gff' 'Bacillus_subtilis.gff')
brev=('EF_chrom' 'EF_pl1' 'EF_pl2' 'EF_pl3' 'LM' 'BS')


# for ((i=0;i<${#termSeqFiles1[@]};++i))
# do
# 	# echo $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles1[i]}.bedgraph
# 	# echo $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles2[i]}.bedgraph
# 	# echo $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles3[i]}.bedgraph
# 	# echo $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${RNAseqFiles[i]}_genomeCoverageBed.bedgraph
# 	# echo $pathToParentDirectory/Termi/Genomes/${gffFiles[i]}
# 	# echo ${brev[i]}

# 	python $pathToParentDirectory/Termi/Scripts/ratioTermRNASeq.py  \
# 			-gff $pathToParentDirectory/Termi/Genomes/${gffFiles[i]} \
# 			-ts $pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles1[i]}.bedgraph \
# 				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles2[i]}.bedgraph \
# 				$pathToParentDirectory/Termi/Alignments/TermSeq/stop_filtered_trimmed_${termSeqFiles3[i]}.bedgraph \
# 			-rs $pathToParentDirectory/Termi/Alignments/RNASeq/sorted_filtered_trimmed_${RNAseqFiles[i]}_genomeCoverageBed.bedgraph \
# 			-bed $pathToParentDirectory/Termi/Alignments/knownTerminators_Bacillus_subtilis.bed \
# 			-c 100 \
# 			-o $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}
# done

# printf "\n##########################################################"
# printf '\n   max. Term-Seq vs. avg. RNA-Seq over all replicates'
# printf "\n##########################################################\n\n"

################################################################################################################################################

# # get fasta files from ratio results (all counts; 120 nucs long), run RNIE, filter out counts with RNIE scores over 20 as "known positives"

fastaFiles=('Enterococcus_faecalis_chromosome.fasta' 'Enterococcus_faecalis_plasmid1.fasta'\
	'Enterococcus_faecalis_plasmid2.fasta' 'Enterococcus_faecalis_plasmid3.fasta'\
	'Listeria_monocytogenes.fasta')
gffFiles2=('Enterococcus_faecalis_Chromosome.gff3' 'Enterococcus_faecalis_plasmid1.gff3'\
	'Enterococcus_faecalis_plasmid2.gff3' 'Enterococcus_faecalis_plasmid3.gff3'\
	'Listeria_monocytogenes.gff')
brev2=('EF_chrom' 'EF_pl1' 'EF_pl2' 'EF_pl3' 'LM')

# for ((i=0;i<${#fastaFiles[@]};++i))
# do
# 	# echo $pathToParentDirectory/Termi/Genomes/${fastaFiles[i]}
# 	echo $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}_50_nucsRNIE_TSvsRS.bed
# 	# echo $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}_50_nucsRNIE_TScountsOverlappingGenes.fasta

# 	bedtools getfasta -fi $pathToParentDirectory/Termi/Genomes/${fastaFiles[i]} \
# 		-bed $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}_50_nucsRNIE_TSvsRS.bed\
#  		-name | \
# 	cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}_50_nucsRNIE_TSvsRS.fasta 
# 	perl $pathToRNIE/bin/rnie.pl \
# 		--gene \
# 		-f $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS.fasta \
# 	 	-md $pathToRNIE/Sequences/ \
# 		-th 0 \
# 	 	-p $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS_th0
# 	python $pathToParentDirectory/Termi/Scripts/filterRNIE.py\
# 		-i $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff \
# 		-o $pathToParentDirectory/Termi/Results/wholeGenome_${brev2[i]}_
# done


# printf "\n##########################################################"
# printf '\n            all RNIE scores filtered'
# printf "\n##########################################################\n\n"

################################################################################################################################################

# # classify counts into "predicted positives/terminators" and "predicted negatives"

# mkdir $pathToParentDirectory"/Termi/Results/Classification"

# # # b.subtilis
# python $pathToParentDirectory/Termi/Scripts/classification.py  \
# 		-pos $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_BS_50_nucs_TScountsOverlappingTerminators \
# 		-neg $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_BS_50_nucs_TScountsOverlappingGenes \
# 		-all $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_BS_50_nucs_TSvsRS \
# 		-gff $pathToParentDirectory/Termi/Genomes/Bacillus_subtilis.gff  \
# 		-term $pathToParentDirectory/Termi/Alignments/knownTerminators_Bacillus_subtilis.bed \
# 		-o $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_

# # # rest
# for ((i=0;i<${#gffFiles2[@]};++i))
# do
# 	python $pathToParentDirectory/Termi/Scripts/classification.py  \
# 		-pos $pathToParentDirectory/Termi/Results/wholeGenome_${brev[i]}_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed \
# 		-neg $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucs_TScountsOverlappingGenes \
# 		-all $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev2[i]}_50_nucs_TSvsRS \
# 		-gff $pathToParentDirectory/Termi/Genomes/${gffFiles2[i]} \
# 		-o $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_
# done

# printf "\n##########################################################"
# printf '\n    classified predicted positives and negatives'
# printf "\n##########################################################\n\n"

################################################################################################################################################

# # removing predictions over 150 nucleotides downstream of genes 

# mkdir $pathToParentDirectory"/Termi/Results/Distance"

# for ((i=0;i<${#gffFiles[@]};++i))
# do
# 	python $pathToParentDirectory/Termi/Scripts/position.py \
# 		-pos $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_${brev[i]}_60_predictedTerminators_NO_knownTerminators_NO_genes.bed \
# 		-neg $pathToParentDirectory/Termi/Results/Classification/wholeGenome_filtered_trim_scaled_${brev[i]}_60_predictedNegatives_NO_knownTerminators_NO_genes.bed \
# 		-gff $pathToParentDirectory/Termi/Genomes/${gffFiles[i]}  \
# 		-o $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_
# done


# printf "\n##########################################################"
# printf '\n       distances to closest genes calculated'
# printf "\n##########################################################\n\n"

# #############################################################################################################################################

# # # removing predicted positives that are too similar to known terminators (BLAST bitscores > 30 )

# # # make blast db of known terminators
makeblastdb -in $pathToParentDirectory/Termi/knownTerminators/true_upper_UtoT_shortIDs.fa \
			-dbtype nucl -parse_seqids -out $pathToBLASTDB/known_terminators_db

mkdir $pathToParentDirectory"/Termi/Results/BLAST"

fastaFiles=('Enterococcus_faecalis_chromosome.fasta' 'Enterococcus_faecalis_plasmid1.fasta'\
	'Enterococcus_faecalis_plasmid2.fasta' 'Enterococcus_faecalis_plasmid3.fasta'\
	'Listeria_monocytogenes.fasta' 'Bacillus_subtilis.fasta')

for ((i=0;i<${#fastaFiles[@]};++i))
do
	# echo $pathToParentDirectory/Termi/Genomes/${fastaFiles[i]}
	# echo $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}_50_nucsRNIE_TSvsRS.bed
	# echo $pathToParentDirectory/Termi/Results/wholeGenome_filtered_trim_scaled_${brev[i]}_50_nucsRNIE_TScountsOverlappingGenes.fasta

	bedtools getfasta \
		-fi $pathToParentDirectory/Termi/Genomes/${fastaFiles[i]} \
		-bed $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed \
 		-name | \
	cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
		then echo $L; read L; echo $L; fi; done > $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta 
	blastn -query $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta \
		-db $pathToBLASTDB/known_terminators_db \
		-word_size 7 \
		-outfmt 6 |
	sort -k12 -rn > $pathToParentDirectory/Termi/Results/BLAST/wholeGenome_filtered_trim_scaled_${brev[i]}_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab
	python $pathToParentDirectory/Termi/Scripts/filterBLAST.py \
		-term  $pathToParentDirectory/Termi/Results/Distance/wholeGenome_filtered_trim_scaled_${brev[i]}_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed \
		-blast $pathToParentDirectory/Termi/Results/BLAST/wholeGenome_filtered_trim_scaled_${brev[i]}_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab \
		-o $pathToParentDirectory/Termi/Results/BLAST/wholeGenome_filtered_trim_scaled_${brev[i]}_
done


printf "\n##########################################################"
printf '\n\t BLAST "predicted" against "known" Terminators'
printf "\n##########################################################\n\n"