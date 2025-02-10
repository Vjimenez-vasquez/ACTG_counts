#!/bin/bash

input=$1
output=$2
grep ">" $input | sed 's/>//g' > input.txt
for line in $(cat input.txt)
do
echo $line > ${line}.seq
seqtk subseq $input ${line}.seq > ${line}.fa
rm ${line}.seq
done 

for b1 in *.fa
do
genome1=$(basename $b1 .fa)
a1=$(grep -v ">" $b1 | grep -E -o "[A]" | wc -l )
a2=$(grep -v ">" $b1 | grep -E -o "[C]" | wc -l )
a3=$(grep -v ">" $b1 | grep -E -o "[T]" | wc -l )
a4=$(grep -v ">" $b1 | grep -E -o "[G]" | wc -l )
a5=$(grep -v ">" $b1 | grep -E -o "[ACTG]" | wc -l )
a6=$(grep -v ">" $b1 | grep -E -o "[ACTGNRYSWKMBDHV-]" | wc -l )
a7=$(grep -v ">" $b1 | grep -E -o "[RYSWKMBDHV]" | wc -l )
a8=$(echo "scale=4; $a5/$a6" | bc )
echo $genome1 $a1 $a2 $a3 $a4 $a5 $a6 $a7 $a8
done > temp.txt
sed '1i sample A C T G ACTG TOTAL ambiguous completeness' temp.txt | tr " " "\t" > actg_counts_${output}.tsv
cat actg_counts_${output}.tsv ;
rm temp.txt ;

cat actg_counts_${output}.tsv | sort -r -n -k8

rm input.txt *.fa 

