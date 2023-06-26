#!/bin/bash

LISTA=("chr15_random" "chr16_random" "chr17_random"  "chr19_random" "chr1_random"  "chr21_random" "chr22_random" "chr2_random" "chr3_random" "chr4_random"  "chr5_h2_hap1" "chr5_random" "chr6_cox_hap1" "chr6_qbl_hap2" "chr6_random" "chr7_random" "chr8_random" "chrX_random")
CAMINHO=/data/AG-BioInf/guestwalther3/human_genome

for CHR in ${LISTA[@]}
do
        echo ${CHR}
	FASTA=${CAMINHO}/chroms_hg18/${CHR}.fa
	if [ -e "${FASTA}" ]; then
		echo "${FASTA} exists. Running code..."
		python3 ${CAMINHO}/extracting/extractor_HumanAS.py -f ${FASTA} -i ${CAMINHO}/KnownAlt/hg18_${CHR}_knownAlt.txt -o ${CAMINHO}/extracting/hg18_${CHR}
	else
		echo "${FASTA} doesnt exist."
	fi
done