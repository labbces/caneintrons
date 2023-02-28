# add argparse para os caminhos dos arquivos.
# mudar os gtf para ter os meso identificadores que o fasta
# python3 -u "/home/bia/sugarcane_introns_local/extracting and trimming/\
# extractor.py" -i /home/bia/sugarcane_introns_local/caneintrons/scripts/\
# trimming_gff_files/trimmed_gtf/SP803280_D_151-161.gtf -f /home/bia/\
# sugarcane_introns/introns_sugarcane/data/Genomes/GCA_002018215.1_CTBE_SP803280_v1.0_genomic.fna\
#  -o output -c /home/bia/sugarcane_introns_local/spliceator_anaysis/\
# contig_size_data.csv

from Bio import SeqIO
from Bio.Seq import Seq
import pyranges as pr
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--gtf_input", type=str, required=True)
parser.add_argument("-f", "--path_to_fasta", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
parser.add_argument("-c", "--contig_table", type=str, required=True,
                    help='table with contig IDs and their lengths')
# parser.add_argument("-x", "--fasta_index", type=str, required=True)
args = parser.parse_args()

# fasta = '/home/bia/sugarcane_introns/old_nov30/data/Genomes/\
# GCA_002018215.1_CTBE_SP803280_v1.0_genomic.fna'
fasta = args.path_to_fasta
# fasta_index = args.fasta_index
inx = args.path_to_fasta+'.genome_inx'
fasta_index = SeqIO.index_db(inx, fasta, 'fasta')

contig_table = args.contig_table
contig_table_pd = pd.read_csv(contig_table, header=None)
contig_lengths = dict(zip(contig_table_pd[0], contig_table_pd[1]))


# Getting gtf file
# gtf = pr.read_gtf('/home/bia/sugarcane_introns/gtf_teste_cana.gtf')
# gtf = pr.read_gtf('/home/bia/sugarcane_introns/old_nov30/\
# 2.TrimmingGenomeFile/SP803280_141-151.gtf')
gtf = pr.read_gtf(args.gtf_input)

# gtf basic operations
# print(gtf)
# print(gtf.describe("transcript_id"))
# print(gtf.features.introns(by="transcript"))

# Introns features pyranges and its basic operations
introns = gtf.features.introns(by="transcript")
# print(introns)
# print(introns.Start)
# print(introns.End)

# Getting only introns that follows: 150 > intron size > 100
introns = (introns[(introns.End - introns.Start) > 99])
introns = (introns[(introns.End - introns.Start) < 151])
print(introns)


filename = args.output

with open(f'{filename}_intron.fa', 'w') as intron_file:
    for seq_id, Start, End, gene_id, strand in \
        zip(introns.df.Chromosome,
            introns.df.Start, introns.df.End,
            introns.df.gene_id, introns.df.Strand):

        seq_full = fasta_index[seq_id].seq

        id_full = f'>{seq_id}.{Start}-{End}-{strand}-{gene_id}'

        start = Start - 200
        end = End + 200

        if start <= 0:
            start = 1
            print(id_full)

        if end >= contig_lengths[seq_id]:
            end = contig_lengths[seq_id]
            print(id_full)

        intron_seq = seq_full[start:end]

        if strand == "+":
            intron_seq = intron_seq
        elif strand == "-":
            seq = Seq(intron_seq)
            intron_seq = seq.reverse_complement()
        else:
            print(
                f'Strand {strand} is not + or - (gene {gene_id} and\
                     seq_id {seq_id}')

        intron_file.write(f'{id_full}\n')
        intron_file.write(f'{intron_seq}\n')
