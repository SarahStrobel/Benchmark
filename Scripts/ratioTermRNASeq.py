# !/bin/bash

pathToParentDirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# samtools view -S -b "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_chromosome.sam" > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_chromosome.bam"
# samtools view -S -b "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid1.sam" > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid1.bam"
# samtools view -S -b "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid2.sam" > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid2.bam"
# samtools view -S -b "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid3.sam" > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid3.bam"


# # #filter out multimapped reads
# samtools view -b -q 70 "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_chromosome.bam"\
# 	 > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_chromosome.bam"

# samtools view -b -q 70 "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid1.bam"\
# 	 > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_plasmid1.bam"

# samtools view -b -q 70 "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid2.bam"\
# 	 > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_plasmid2.bam"

# samtools view -b -q 70 "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid3.bam"\
# 	 > "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_plasmid3.bam"


# awk '/^@/ || $5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_chromosome.sam" \
# 	> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_chromosome.sam"

# awk '/^@/ || $5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid1.sam" \
# 	> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_plasmid1.sam"

# awk '/^@/ || $5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid2.sam" \
# 	> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_plasmid2.sam"

# awk '/^@/ || $5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/trim_ERR1248404_plasmid3.sam" \
# 	> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_plasmid3.sam"


# samtools sort "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_chromosome.bam" \
# 				-o "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_chromosome.bam" \
# 				> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_chromosome.bam" 

# samtools sort "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_plasmid1.bam" \
# 				-o "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid1.bam" \
# 				> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid1.bam" 

# samtools sort "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_plasmid2.bam" \
# 				-o "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid2.bam" \
# 				> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid2.bam" 

# samtools sort "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_ERR1248404_plasmid3.bam" \
# 				-o "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid3.bam" \
# 				> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid3.bam" 


# samtools index "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_chromosome.bam"
# samtools index "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid1.bam" 
# samtools index "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid2.bam" 
# samtools index "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid3.bam" 



# genomeCoverageBed -d -split -ibam "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_chromosome.bam" \
# 					> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_chromosome_genomeCoverageBed.bedgraph"

# genomeCoverageBed -d -split -ibam "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid1.bam" \
# 					> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid1_genomeCoverageBed.bedgraph"

# genomeCoverageBed -d -split -ibam "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid2.bam" \
# 					> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid2_genomeCoverageBed.bedgraph"

# genomeCoverageBed -d -split -ibam "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid3.bam" \
# 					> "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid3_genomeCoverageBed.bedgraph"


# wc -l "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_chromosome_genomeCoverageBed.bedgraph"
# wc -l "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid1_genomeCoverageBed.bedgraph"
# wc -l "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid2_genomeCoverageBed.bedgraph"
# wc -l "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid3_genomeCoverageBed.bedgraph"


# # # calculate unique reads and append header to BEDGRAPH file
# uniqueReads="$(awk '$5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_chromosome.sam" | wc -l)"
# sed -i '1itrack type=bedGraph name=\"'$uniqueReads'\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20' \
# 	"/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_chromosome_genomeCoverageBed.bedgraph"

# uniqueReads="$(awk '$5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_plasmid1.sam" | wc -l)"
# sed -i '1itrack type=bedGraph name=\"'$uniqueReads'\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20' \
# 	"/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid1_genomeCoverageBed.bedgraph"

# uniqueReads="$(awk '$5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_plasmid2.sam" | wc -l)"
# sed -i '1itrack type=bedGraph name=\"'$uniqueReads'\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20' \
# 	"/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid2_genomeCoverageBed.bedgraph"

# uniqueReads="$(awk '$5 > "0"' "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtered_trim_ERR1248404_plasmid3.sam" | wc -l)"
# sed -i '1itrack type=bedGraph name=\"'$uniqueReads'\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20' \
# 	"/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/filtererd_trim_sorted_ERR1248404_plasmid3_genomeCoverageBed.bedgraph"




# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/EndnucCount/5primePosfromSam_2_Bedgraph_withoutSoft_10102019.py"\
# 										-termSeq "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/" \
# 										-o "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10172019.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248401_chromosome.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248402_chromosome.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248403_chromosome.bedgraph" \
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/chromosome_and_plasmids/filtererd_trim_sorted_ERR1248404_chromosome_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-g 1 500000\
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10172019.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248401_chromosome.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248402_chromosome.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248403_chromosome.bedgraph" \
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/chromosome_and_plasmids/filtererd_trim_sorted_ERR1248404_chromosome_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10172019.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF2.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248401_plasmid2.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248402_plasmid2.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/EF/chromosome_and_plasmids/filtered/stop_ERR1248403_plasmid2.bedgraph" \
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/chromosome_and_plasmids/filtererd_trim_sorted_ERR1248404_plasmid2_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled"




# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/E_faecalis_chromosome.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 

# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/E_faecalis_chromosome.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 




# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid3.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 



# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"

# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"

# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"

# #filter rnie output (everything over bitscore 20 --> positive)
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_"





# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_09232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3"  \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/first500k_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_09232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3"  \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10172019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF2.gff3"  \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled"




# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/first500k_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/first500k_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3"  \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3"  \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF3.gff3"  \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_"



# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/E_faecalis_chromosome.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		

# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/E_faecalis_chromosome.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		


# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid3.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"

# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"



# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/"


# sh ~/Masterarbeit/Skripte/reverseComplement/reverseComplement_10212019.sh


# RNIE plots
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/RNIE_first500k_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/RNIE_first500k_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/first500k_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/first500k_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/first500k_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/"



# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/first500k_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/first500k_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/first500k_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/first500k_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/EF/Chromosome/first500k_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/EF/Chromosome/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/EF/Plasmid3/wholeGenome_"




#######################################################################################################################################3

#bsub
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_2.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERX1304415.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERX1320300.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERX1320301.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/BS/filtererd_trim_sorted_ERX1320302_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/BS/wholeGenome_filtered_trim_scaled"

# #lmon 
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_2.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/GCF_000196035.1_ASM19603v1_genomic_Lmon.gff" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248436.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248437.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248438.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/LM/filtered_trim_sorted_ERR1248439_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled"

# # # #spneu
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_SP.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Streptococcus_pneumoniae_tigr4_ASM688v1_37_new.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_SRR7160964.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_SRR7160965.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_SRR7160966.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_SRR7160967.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/SP/combined_316_RNASeqND0min_genomeCoverage.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled"

# #efaec
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_2.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248401_chromosome.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248402_chromosome.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248403_chromosome.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/chromosome_and_plasmids/filtererd_trim_sorted_ERR1248404_chromosome_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_2.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF1.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248401_plasmid1.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248402_plasmid1.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248403_plasmid1.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/chromosome_and_plasmids/filtererd_trim_sorted_ERR1248404_plasmid1_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled"



# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_2.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF2.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248401_plasmid2.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248402_plasmid2.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248403_plasmid2.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/chromosome_and_plasmids/filtererd_trim_sorted_ERR1248404_plasmid2_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_2.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF3.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248401_plasmid3.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248402_plasmid3.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_filtered_trim_ERR1248403_plasmid3.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/EF/chromosome_and_plasmids/filtererd_trim_sorted_ERR1248404_plasmid3_genomeCoverageBed.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled"



###########################################################################################################
# #lmon 
# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/GCF_000196035.1_ASM19603v1_genomic_Lmon.fna" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 

# #spneu
# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/NC_003028.3_S_pneu.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 


# #efaec
# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/E_faecalis_chromosome.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 

# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid1.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 

# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid2.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 

# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid3.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 


# #lmon 
# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" \
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"


# #spneu
# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th10"


# #efaec
# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"


# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"

# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"

# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0"



#filter rnie output (everything over bitscore 20 --> positive)

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th10-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th0-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_"


###########################################################################################################

# #bsub
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/BS/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingTerminators" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/BS/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/BS/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene  "/scr/k70san3/stsarah/Genomes/GCF_000009045.1_ASM904v1_genomic.gff" \
# 						-term $pathToParentDirectory"/Alignments/knownTerminators.bed"\
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/BS/wholeGenome_filtered_trim_scaled"

# #lmon
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/GCF_000196035.1_ASM19603v1_genomic_Lmon.gff"  \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/LM/wholeGenome_filtered_trim_scaled"

# #spneu
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Streptococcus_pneumoniae_tigr4_ASM688v1_37_new.gff3"  \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/SP/wholeGenome_filtered_trim_scaled"

# #efaec
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3"  \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF1.gff3"   \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF2.gff3"   \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF3.gff3"   \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled"


#####################################################################################################################################

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/BS/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/BS/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/GCF_000009045.1_ASM904v1_genomic.gff"  \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/LM/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/LM/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/GCF_000196035.1_ASM19603v1_genomic_Lmon.gff"   \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/SP/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/SP/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Streptococcus_pneumoniae_tigr4_ASM688v1_37_new.gff3"   \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_Chromosome.gff3"   \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF1.gff3"    \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF2.gff3"    \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_"

# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Enterococcus_faecalis_v583_ASM778v1_45_chromosome_pTEF3.gff3"   \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_"

##########################################################################################################################################
# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/Bacillus_subtilis.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		

# #lmon
# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/GCF_000196035.1_ASM19603v1_genomic_Lmon.fna" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 	

# #spneu
# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/NC_003028.3_S_pneu.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		

# #efaec
# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/E_faecalis_chromosome.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 	


# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid1.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		


# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid2.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 	


# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/E_faecalis_plasmid3.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		

# #############################################################################
# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"



# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_"

# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_"

# #############################################################################################

# # sh ~/Masterarbeit/Skripte/reverseComplement/reverseComplement_10232019.sh


# # RNIE plots
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/BS/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/BS/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/BS/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/BS/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/BS/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/BS/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/LM/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/LM/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/LM/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/LM/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/LM/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/SP/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/SP/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/SP/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/SP/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/SP/wholeGenome_"




# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid1/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid1/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid2/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid2/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_09232019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_"

# ################################################################################
# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/BS/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/BS/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/BS/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/BS/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/BS/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/LM/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/LM/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/LM/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/LM/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/LM/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/SP/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/SP/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/SP/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/SP/wholeGenome_"




# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Chromosome/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/EF/Chromosome/wholeGenome_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid1/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/EF/Plasmid1/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid2/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/EF/Plasmid2/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10212019/EF/Plasmid3/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10212019/EF/Plasmid3/wholeGenome_"


###################################################################################################


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/EndnucCount/5primePosfromSam_2_Bedgraph_withoutSoft_10242019.py"\
# 		-rnaSeq "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/SP/combined_316_RNASeqND0min_genomeCoverage.bedgraph"\
# 		-termSeq "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/" \
# 		-o "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/stop_new_"

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/EndnucCount/5primePosfromSam_2_Bedgraph_withoutSoft_10102019.py"\
# 		-termSeq "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/" \
# 		-o "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/stop_new_"


# genomeCoverageBed -bga -split -ibam "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/filtered_trim_sorted_SRR7160967.bam" \
# 					> "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/filtered_trim_sorted_SRR7160967_genomeCoverageBed.bedgraph"

###################################################################################################

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Ratio/ratioTermRNASeq_10222019_SP.py"  \
# 						-gff "/scr/k70san3/stsarah/Genomes/Streptococcus_pneumoniae_tigr4_genbank2gff3.gff3" \
# 						-ts "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_new_filtered_trim_SRR7160964_.bedgraph" \
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_new_filtered_trim_SRR7160965_.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_new_filtered_trim_SRR7160966_.bedgraph"\
# 							"/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq_all/stop_new_filtered_trim_SRR7160967_.bedgraph"\
# 						-rs "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/SP/combined_316_RNASeqND0min_genomeCoverage.bedgraph" \
# 						-bed $pathToParentDirectory"/Alignments/knownTerminators.bed" \
# 						-c 100 \
# 						-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled"

# bedtools getfasta -fi "/scr/k70san3/stsarah/Genomes/NC_003028.3_S_pneu.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.bed"\
#  		-name | \
# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta" 


# perl "/homes/brauerei/stsarah/Masterarbeit/Infernal/bin/rnie.pl"\
# 		 --gene -f "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS.fasta"\
# 		 -md "/homes/brauerei/stsarah/Masterarbeit/Infernal/Sequences/" \
# 		 -th 0 \
# 		 -p "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th10"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/filterRNIE.py" \
# 	-i "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucsRNIE_TSvsRS_th10-geneMode-rnie.gff" \
# 	-o "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_"


# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/Classification/SVM_gradStudAlgo_combined_10232019.py"  \
# 						-pos "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_50_nucsRNIE_TSvsRS_RNIEover20.bed" \
# 						-neg "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TScountsOverlappingGenes" \
# 						-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 						-gene "/scr/k70san3/stsarah/Genomes/Streptococcus_pneumoniae_tigr4_genbank2gff3.gff3"  \
# 						-o "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10242019/SP/wholeGenome_filtered_trim_scaled"


# python $pathToParentDirectory"/Scripts/position_everything.py" \
# 		-pos "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10242019/SP/wholeGenome_filtered_trim_scaled60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-neg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10242019/SP/wholeGenome_filtered_trim_scaled60_predictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-gene "/scr/k70san3/stsarah/Genomes/Streptococcus_pneumoniae_tigr4_ASM688v1_37_new.gff3"   \
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_"

# bedtools getfasta -fi  "/scr/k70san3/stsarah/Genomes/NC_003028.3_S_pneu.fasta" \
# 		-bed "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
#  		-name | \

# cat | while read L; do if [[ $L == *+ ]]; then echo $L; read L; echo $L | rev | tr "ATGC" "TACG" ; elif [[ $L == *- ]]; \
# 		then echo $L; read L; echo $L; fi; done > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" 		


# blastn -query "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.fasta" \
# 		-db ~/blastdb/new_known_terminators_db \
# 		-word_size 7 \
# 		-outfmt 6 |
# sort -k12 -rn > "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"


# python $pathToParentDirectory"/Scripts/filterBLAST_09232019.py" \
# 		-term  "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_Distance_60_predictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-blast "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_BLAST_60_predictedTerminators_NO_knownTerminators_NO_genes.tab"\
# 		-o "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_"

# sh ~/Masterarbeit/Skripte/reverseComplement/reverseComplement_10242019.sh

# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_scores/RNIE_SVM_gradStudAlgo_10252019.py"\
# 		-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS"\
# 		-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10242019/SP/RNIE_wholeGenome_positives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10242019/SP/RNIE_wholeGenome_negatives_NO_knownTerminators_NO_genes_th0-geneMode-rnie.gff" \
# 		-bedPos "//scr/k70san3/stsarah/Output_sklearnSVM_500k/10242019/SP/wholeGenome_filtered_trim_scaledpredictedTerminators_NO_knownTerminators_NO_genes.bed"\
# 		-bedNeg "/scr/k70san3/stsarah/Output_sklearnSVM_500k/10242019/SP/wholeGenome_filtered_trim_scaledpredictedNegatives_NO_knownTerminators_NO_genes.bed"\
# 		-o "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10242019/SP/wholeGenome_"



# python "/homes/brauerei/stsarah/Masterarbeit/Skripte/RNIE_Position/RNIE_Position_combined_09242019.py" \
# 	-all "/scr/k70san3/stsarah/Output_ratioTermRNASeq/10242019/SP/wholeGenome_filtered_trim_scaled_50_nucs_TSvsRS" \
# 	-rniePos "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10242019/SP/wholeGenome_RNIE_scores_predictedPositves_NO_knownTerminators_NO_genes.bed" \
# 	-rnieNeg "/scr/k70san3/stsarah/Output_RNIE_SVM_gradStudAlgo/10242019/SP/wholeGenome_RNIE_scores_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distNeg "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_Distance_predictedNegatives_NO_knownTerminators_NO_genes.bed" \
# 	-distPos "/scr/k70san3/stsarah/Output_Position_gradStudAlgo/10242019/SP/wholeGenome_filtered_trim_scaled_Distance_predictedTerminators_NO_knownTerminators_NO_genes.bed" \
# 	-o "/scr/k70san3/stsarah/Output_RNIE_Position_combined/10242019/SP/wholeGenome_"



python "/homes/brauerei/stsarah/Masterarbeit/Skripte/EndnucCount/5primePosfromSam_2_Bedgraph_withoutSoft_10252019.py"\
		-rnaSeq "/scr/k70san3/stsarah/Novoalign_Alignments/RNA_Seq/SP/combined_316_RNASeqND0min_genomeCoverage.bedgraph"\
		-termSeq "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/" \
		-o "/scr/k70san3/stsarah/Novoalign_Alignments/Term_Seq/SP/filtered/stop_10252019_"