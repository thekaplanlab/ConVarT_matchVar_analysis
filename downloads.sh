#!/bin/bash

# Download latest wormbase variant file

while read line2; do
    while IFS= read -r line; do
        if [[ $line2 != $line ]]; then
            curl -O ftp://ftp.wormbase.org/pub/wormbase/releases/${line}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.${line}.annotations.gff3.gz
            zcat c_elegans.PRJNA13758.${line}.annotations.gff3.gz | grep variation=WBVar | grep aachange= > wormbase_variants_aachange.txt
            curl -o c_elegans_protein_wb_fa.gz ftp://ftp.wormbase.org/pub/wormbase/releases/${line}/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.${line}.protein.fa.gz
            gunzip c_elegans_protein_wb.fa.gz
            echo "$line"> wb_version2
        fi
    done < wb_version
done < wb_version2
