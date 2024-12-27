# ACTG_counts
actg_counts

## Just run for every file.fasta 
```r
#!/bin/bash

grep ">" file.fasta | sed 's/>//g' > header.txt ;

var=$(seq 2 2 16768)
for i in ${var[@]}
do
sed -n "${i}p" file.fasta | grep -E -o "[ACTG]" | wc -l
done > actg.txt

paste -d ',' header.txt actg.txt | sed '1i seq,actg '> counts.txt

rm header.txt actg.txt ;
cat counts.txt ; 

echo "############ " ; 
echo "## DONE ! ## " ;
echo "############ " ;

exit
```

## Just run for every file.fas (improved)
```r
#!/bin/bash

for b1 in *.fa
do
genome1=$(basename $b1 .fa)
a1=$(grep -v ">" $b1 | grep -E -o "[A]" | wc -l )
a2=$(grep -v ">" $b1 | grep -E -o "[C]" | wc -l )
a3=$(grep -v ">" $b1 | grep -E -o "[T]" | wc -l )
a4=$(grep -v ">" $b1 | grep -E -o "[G]" | wc -l )
a5=$(grep -v ">" $b1 | grep -E -o "[ACTG]" | wc -l )
a6=$(grep -v ">" $b1 | grep -E -o "[ACTGN-]" | wc -l )
a7=$(echo "scale=4; $a5/$a6" | bc )
echo $genome1 $a1 $a2 $a3 $a4 $a5 $a6 $a7
done > temp.txt
sed '1i sample A C T G ACTG TOTAL completeness' temp.txt | tr " " "\t" > actg_counts.tsv
cat actg_counts.tsv ;
rm temp.txt ;

cat actg_counts.tsv | sort -r -n -k8
```

## the last code
```r
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
a1=$(grep -v ">" $b1 | grep -E -o "[a]" | wc -l )
a2=$(grep -v ">" $b1 | grep -E -o "[c]" | wc -l )
a3=$(grep -v ">" $b1 | grep -E -o "[t]" | wc -l )
a4=$(grep -v ">" $b1 | grep -E -o "[g]" | wc -l )
a5=$(grep -v ">" $b1 | grep -E -o "[actg]" | wc -l )
a6=$(grep -v ">" $b1 | grep -E -o "[actgn-]" | wc -l )
a7=$(echo "scale=4; $a5/$a6" | bc )
echo $genome1 $a1 $a2 $a3 $a4 $a5 $a6 $a7
done > temp.txt
sed '1i sample A C T G ACTG TOTAL completeness' temp.txt | tr " " "\t" > actg_counts_${output}.tsv
cat actg_counts_${output}.tsv ;
rm temp.txt ;

cat actg_counts_${output}.tsv | sort -r -n -k8

rm input.txt *.fa 
```

## usage: 
```r
./actg_counts_improved input.fasta output_prefix
```

## output: 
```r
actg_counts_output_prefix.tsv
```
