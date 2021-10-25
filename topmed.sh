#!/bin/bash

# Script to download and process TOPMed variants.

mkdir topmed
for i in {1..22}; do wget http://ftp.ensembl.org/pub/grch37/release-104/variation/vcf/homo_sapiens/homo_sapiens-chr${i}.vcf.gz -P topmed/ & done

wget http://ftp.ensembl.org/pub/grch37/release-104/variation/vcf/homo_sapiens/homo_sapiens-chrX.vcf.gz -P topmed/
wget http://ftp.ensembl.org/pub/grch37/release-104/variation/vcf/homo_sapiens/homo_sapiens-chrY.vcf.gz -P topmed/

for i in ~/topmed/*.vcf.gz; do date; echo $i; ./vep --cache --offline --fork 30 --buffer_size 50000 --force_overwrite --no_stats --assembly GRCh37 --protein --symbol -sift b -polyphen b --gene_phenotype -i $i -o "${i%.vcf.gz}_1.txt"; done

awk 'NR==47' homo_sapiens-chr4_1.txt > all_output.txt
awk 'FNR>47 && ($10 != "-") {print}' *_1.txt >> all_output.txt
awk -F"SIFT=" '{print $2}' all_output.txt | cut -f1 -d";" > sift.txt
awk -F"PolyPhen=" '{print $2}' all_output.txt | cut -f1 -d";" > polyphen.txt
awk -F"ENSP=" '{print $2}' all_output.txt | cut -f1 -d";" > ensp.txt
paste all_output.txt ensp.txt sift.txt polyphen.txt > all_output_clean.txt
awk '{print $1, $4, $5, $7, $10, $11, $15, $16, $17}' all_output_clean.txt > all_output_clean1.txt

