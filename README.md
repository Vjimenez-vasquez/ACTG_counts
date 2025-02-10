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

## Now usea FILTER1 to extract the contigs names over a given length (50000, 100000 or whatever)
```r
filter1 <- function(extension,prefix2,title,contig.length){
  
  extension <- extension
  prefix2 <- prefix2
  title <- title
  contig.length <- as.numeric(contig.length)
  
  r <- grep(extension,dir(),value=T) ;
  
  r1 <- data.frame() ; r2 <- 0 ; r3 <- 0 ; r4 <- 0 ; r5 <- 0 ; r6 <- 0 ; 
  r7 <- 0 ; r8 <- 0 ; r9 <- 0 ; r10 <- 0 ; r11 <- 0 ; r12 <- 0 ; r13 <- 0 ;
  r14 <- 0 ; r15 <- 0 ;
  
  for (i in 1:length(r)){
    dim(r1)
    r1 <- read.csv(r[i], header=T, sep="\t")
    r2 <- gsub(extension,"",r[i],fixed=F)
    r3 <- rep(gsub(prefix2,"",r2,fixed=F),nrow(r1))
    r4 <- append(r4,r3)
    r5 <- append(r5,r1$sample)
    r6 <- append(r6,r1$A)
    r7 <- append(r7,r1$C)
    r8 <- append(r8,r1$T)
    r9 <- append(r9,r1$G)
    r10 <- append(r10,r1$ACTG)
    r11 <- append(r11,r1$TOTAL)
    r12 <- append(r12,r1$ambiguous)
    r13 <- append(r13,r1$completeness)
    
    r14 <- paste0(r3,"@",rep("contig",nrow(r1)),1:nrow(r1),"_",r1$TOTAL,"pb")
    r15 <- append(r15, r14)
  }
  
  #check the items# 
  s1 <- data.frame(r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r15)
  s2 <- s1[2:nrow(s1),]
  names(s2) <- c("sample","original.contig","A","C","T","G","ACTG","TOTAL","ambiguous","completeness","suggested.contig.name")
  
  write.table(s2,paste0(title,".total",".tsv"), row.names=F, col.names=T, quot=F, sep=" \t")
  
  s3 <- s2[s2$TOTAL > contig.length, ]
  write.table(s3,paste0(title,".",as.character(contig.length),".tsv"), row.names=F, col.names=T, quot=F, sep=" \t")
  
  t1 <- 0 ; 
  for(j in unique(s3$sample)){
    t1 <- s3[s3$sample %in% j, 2]
    write.table(t1, paste0(j,".",as.character(contig.length),".contigs.txt"), row.names=F, col.names=F, quot=F)
  }
  
}

## EXAMPLE ##
## filter1(".tsv","actg_counts_","pilon","1000000") ##

## extension: extension de los archivos a procesar (".txt",".tsv"), es tambie el ultimo patron a remover para obtener el prefijo de las muestras 
## prefix2: palabar inicial a remover para obtener el prefijo de las muestras
## title: titulo de los archivos de salida
## contig.length: longitud minima de contig 
```
