#!/bin/bash
LISTA=("chr1" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr15_random" "chr16" "chr16_random" "chr17" "chr17_random" "chr18" "chr19" "chr19_random" "chr1_random" "chr2" "chr20" "chr21" "chr21_random" "chr22" "chr22_random" "chr2_random" "chr3" "chr3_random" "chr4" "chr4_random" "chr5" "chr5_h2_hap1" "chr5_random" "chr6" "chr6_cox_hap1" "chr6_qbl_hap2" "chr6_random" "chr7" "chr7_random" "chr8" "chr8_random" "chr9" "chrX" "chrX_random" "chrY")

for CHR in "${LISTA[@]}"
do
    echo "Starting ${CHR}"
    python3 -u "/home/bia/sugarcane_introns_local/known_alt_separator.py" -c ${CHR} >> hg18_${CHR}_knownAlt.txt
done