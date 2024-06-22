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