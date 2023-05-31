#!/bin/bash

GTF_FILES='/home/bia/sugarcane_introns_local/caneintrons/scripts/trimming_gff_files/trimmed_gtf/SP803280_D_*'
GENOME_FASTA='/home/bia/sugarcane_introns/introns_sugarcane/data/Genomes/GCA_002018215.1_CTBE_SP803280_v1.0_genomic.fna'
CONTIG_TABLE='/home/bia/sugarcane_introns_local/spliceator_anaysis/contig_size_data.csv'

for GTF_FILE in $GTF_FILES
do
    IN=`basename $GTF_FILE`
    OUTPUT=${IN/.GTF/200}
    python3 /home/bia/sugarcane_introns_local/caneintrons/scripts/extracting_introns/extractor.py -i ${GTF_FILE} -f ${GENOME_FASTA} -o ${OUTPUT} -c ${CONTIG_TABLE} 
done

