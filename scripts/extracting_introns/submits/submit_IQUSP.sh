#!/bin/bash

GTF_FILES='/home/bia/sugarcane_introns_local/extracting_and_trimming/IQ_USP_trimmed_gtf/Souza_*.gtf'
GENOME_FASTA='/home/bia/sugarcane_introns_local/data/Genomes/Souza.sc.mlc.cns.sgl.utg.cga7.importdb_id_changed.fna'
CONTIG_TABLE='/home/bia/sugarcane_introns_local/extracting_and_trimming/Souza_contig_table.csv'

for GTF_FILE in $GTF_FILES
do
    IN=`basename $GTF_FILE`
    OUTPUT=${IN/.gtf/200}
    python3 /home/bia/sugarcane_introns_local/caneintrons/scripts/extracting_introns/extractor.py -i ${GTF_FILE} -f ${GENOME_FASTA} -o ${OUTPUT} -c ${CONTIG_TABLE} >> IQUSP_extrator_output.txt
done