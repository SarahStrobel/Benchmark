# !/bin/bash


pathToParentDirectory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

mkdir $pathToParentDirectory/"Benchmark"

brev=('BS' 'EF_chrom' 'LM' 'SP' )

mkdir $pathToParentDirectory/Benchmark/RNIE/
mkdir $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive
mkdir $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive
mkdir $pathToParentDirectory/Benchmark/RNAmotif/
mkdir $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive
mkdir $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive
mkdir $pathToParentDirectory/Benchmark/iTerm_PseKNC/
mkdir $pathToParentDirectory/Benchmark/iTerm_PseKNC/1NegativePerPositive

for file in $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/*.fasta
do
	echo $file
	if /usr/bin/time rnie.pl \
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

for file in $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/*.fasta
do
	echo $file
	if /usr/bin/time rnie.pl \
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

cd $pathToParentDirectory/Termi/iTerm-PseKNC_modified/

for file in $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/*.csv
do
	echo $file
	if /usr/bin/time python $pathToParentDirectory/Termi/iTerm-PseKNC_modified/iTerm-PseKNC_modified.py \
		$file \
		$pathToParentDirectory/Benchmark/iTerm_PseKNC/1NegativePerPositive/1_predicted_"$(basename "$file" .csv)"; then
			echo 'worked'
	else
		exit 1
	fi
done

######################################################################
######################################################################



for file in $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/*.gff
do
	echo $file
	perl $pathToParentDirectory/Termi/RNIE/rnie2gff.pl -r $file
done


for file in $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/*.out
do
	perl $pathToParentDirectory/Termi/RNIE/terminator-lesnik2gff.pl -t $file 
done

or file in $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive/*.gff
do
	echo $file
	perl $pathToParentDirectory/Termi/RNIE/rnie2gff.pl -r $file
done


for file in $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive/*.out
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
	python $pathToParentDirectory/Termi/Scripts/iTerm2gff.py \
	-iterm ${resultFiles[i]} \
	-i ${fastaFiles[i]} \
	-o $pathToParentDirectory/Benchmark/iTerm_PseKNC/iTerm_results_from_wl/"$(basename "${resultFiles[i]}" .txt)"
done

#######################################################################
#######################################################################

mkdir $pathToParentDirectory/"Benchmark/ROC"
mkdir $pathToParentDirectory/"Benchmark/ROC/1NegativePerPositive"
mkdir $pathToParentDirectory/"Benchmark/ROC/100NegativesPerPositive"

for ((i=0;i<${#brev[@]};++i))
do
	a=$(cat $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)
	b=$(cat $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)

	python $pathToParentDirectory/Termi/Scripts/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/iTerm_PseKNC/iTerm_results_from_wl/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled_result.gff \
	-false $pathToParentDirectory/Benchmark/iTerm_PseKNC/iTerm_results_from_wl/1_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled_result.gff \
	-o $pathToParentDirectory/Benchmark/ROC/1NegativePerPositive/iTerm_${brev[i]}_shuffled_ \
	-nucs $((a+b))

	python $pathToParentDirectory/Termi/Scripts/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled_th0-geneMode-rnie.bits.gff \
	-false $pathToParentDirectory/Benchmark/RNIE/1NegativePerPositive/1_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled_th0-geneMode-rnie.bits.gff  \
	-o $pathToParentDirectory/Benchmark/ROC/1NegativePerPositive/RNIEth0_${brev[i]}_shuffled_ \
	-nucs $((a+b))

	python $pathToParentDirectory/Termi/Scripts/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/1NegativePerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/1_predicted_${brev[i]}_embedded_predictedTerminators_shuffled.terminator-lesnik.out.dG_score.gff\
	-false $pathToParentDirectory/Benchmark/RNAmotif/1NegativePerPositive/1_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled.terminator-lesnik.out.dG_score.gff \
	-o $pathToParentDirectory/Benchmark/ROC/1NegativePerPositive/RNAmotif_${brev[i]}_shuffled_ \
	-nucs $((a+b))
done

for ((i=0;i<${#brev[@]};++i))
do
	a=$(cat $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)
	b=$(cat $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta | grep -v \> | tr -d "\r\n"|wc -c)

	python $pathToParentDirectory/Termi/Scripts/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive/100_predicted_${brev[i]}_embedded_predictedTerminators_shuffled_th0-geneMode-rnie.bits.gff \
	-false $pathToParentDirectory/Benchmark/RNIE/100NegativesPerPositive/100_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled_th0-geneMode-rnie.bits.gff  \
	-o $pathToParentDirectory/Benchmark/ROC/100NegativesPerPositive/RNIEth0_${brev[i]}_shuffled_ \
	-nucs $((a+b))

	python $pathToParentDirectory/Termi/Scripts/computeROC.py \
	-pos $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_predictedTerminators_shuffled.fasta \
	-neg $pathToParentDirectory/Termi/Results/Embedded/100NegativesPerPositive/${brev[i]}_embedded_shuffledNegatives_shuffled.fasta \
	-true $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive/100_predicted_${brev[i]}_embedded_predictedTerminators_shuffled.terminator-lesnik.out.dG_score.gff\
	-false $pathToParentDirectory/Benchmark/RNAmotif/100NegativesPerPositive/100_predicted_${brev[i]}_embedded_shuffledNegatives_shuffled.terminator-lesnik.out.dG_score.gff \
	-o $pathToParentDirectory/Benchmark/ROC/100NegativesPerPositive/RNAmotif_${brev[i]}_shuffled_ \
	-nucs $((a+b))
done

#######################################################################
#######################################################################

