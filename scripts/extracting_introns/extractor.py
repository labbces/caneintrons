# add argparse para os caminhos dos arquivos.
# mudar os gtf para ter os meso identificadores que o fasta

from Bio import SeqIO
import pyranges as pr
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--gtf_input", type=str, required=True)
parser.add_argument("-f", "--path_to_fasta", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
#parser.add_argument("-x", "--fasta_index", type=str, required=True)
args = parser.parse_args()

#fasta = '/home/bia/sugarcane_introns/old_nov30/data/Genomes/GCA_002018215.1_CTBE_SP803280_v1.0_genomic.fna'
fasta = args.path_to_fasta
#fasta_index = args.fasta_index
inx = args.path_to_fasta+'.genome_inx'
fasta_index = SeqIO.index_db(inx, fasta, 'fasta')


# Getting gtf file
# gtf = pr.read_gtf('/home/bia/sugarcane_introns/gtf_teste_cana.gtf')
# gtf = pr.read_gtf('/home/bia/sugarcane_introns/old_nov30/2.TrimmingGenomeFile/SP803280_141-151.gtf')
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
    for seq_id, Start, End, gene_id in zip(introns.df.Chromosome, introns.df.Start, introns.df.End, introns.df.gene_id):
        seq_full = fasta_index[seq_id].seq

        id_full = f'>{seq_id}.{Start}-{End}-{gene_id}'
        start = Start - 200
        end = End + 200
        intron_seq = seq_full[start:end]

        intron_file.write(f'{id_full}\n')
        intron_file.write(f'{intron_seq}\n')
