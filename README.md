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
```


