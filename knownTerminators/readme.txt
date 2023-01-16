The file original-true.fa is a FASTA file of known terminators from the RNIE software package.  Its source URL is https://raw.githubusercontent.com/ppgardne/RNIE-benchmark/master/training-test/true.fa

The file true.fa is the same thing, but the sequences were made upper case and U->T.  Like so:

cat original-true.fa | sed "s/[Uu]/t/g" | awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1 " " $2}}' > true.fa

The file dependent_data2.fa is from the Term-PseKNC software.  Its source URL is http://lin-group.cn/server/iTerm-PseKNC/dependent_data2.csv
