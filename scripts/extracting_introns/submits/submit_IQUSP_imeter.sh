#!/bin/bash

GTF_FILES='/home/bia/sugarcane_introns_local/extracting_and_trimming/IQ_USP_trimmed_gtf/Souza_*.gtf'
GENOME_FASTA='/home/bia/sugarcane_introns_local/data/Genomas/Souza.sc.mlc.cns.sgl.utg.cga7.importdb_id_changed.fna'
CONTIG_TABLE='/home/bia/sugarcane_introns_local/data/Contig_and_Strand_tables/Souza_contig_table.csv'

for GTF_FILE in $GTF_FILES
do
    IN=`basename $GTF_FILE`
    OUTPUT=${IN/.gtf/"_imeter"}
    python3 /home/bia/sugarcane_introns_local/extractor_V2.py -i ${GTF_FILE} -f ${GENOME_FASTA} -o ${OUTPUT} -c ${CONTIG_TABLE} -t IMEter
done






