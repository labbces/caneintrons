# script applied to extract intron + flanking region sequences from genome file
# and gtf file using pyranges library
# gff files need to be converted into gft using for exemple gtf read
# It extracts introns/splice sites flanking reagions based on IMEter, spliceator and splice2deep requirements
# intron size is between 100-150bp

from Bio import SeqIO
from Bio.Seq import Seq
import pyranges as pr
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--gtf_input", type=str, required=True)
parser.add_argument("-f", "--path_to_fasta", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True,
                    help='for -o XXXX, output file would be XXXX_[donor/acceptor]_intron.fa')
parser.add_argument("-c", "--contig_table", type=str, required=True,
                    help='table with contig IDs and their lengths - generated with get_contig_table.py')
parser.add_argument("-t", "--trimming_type", choices=['IMEter', 'spliceator',
                    'splice2deep'], help='Options: IMEter SSXXXXXXSS, spliceator 200XSSXXXXXXX200X,\
                         splice2deep for donor and acceptor 300SXS300X. SS = splice site X = nucleotides\
                        ', default='IMEter', required=True)
parser.add_argument("-b", "--flanking_size", type=int,
                    help='Set specific flanking regions. Defaults: IMEter = 0, Spliceator = 200, Splice2Deep = 300', default=None, required=False)
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

filename = args.output
filename_donor = filename + "_donor"
filename_acceptor = filename + "_acceptor"

gtf = pr.read_gtf(args.gtf_input)

flanking_size_dict = {'IMEter': 0, 'spliceator': 200, 'splice2deep': 300}

# Introns features pyranges and its basic operations
introns = gtf.features.introns(by="transcript")
# print(introns)
# print(introns.Start)
# print(introns.End)

# Getting only introns that follows: 150 > intron size > 100
introns = (introns[(introns.End - introns.Start) > 99])
introns = (introns[(introns.End - introns.Start) < 151])
# print(introns)


trimming_type = args.trimming_type

with open(f'{filename}_intron.fa', 'w') as intron_file:
    for seq_id, Start, End, gene_id, strand in \
        zip(introns.df.Chromosome,
            introns.df.Start, introns.df.End,
            introns.df.gene_id, introns.df.Strand):

        seq_full = fasta_index[seq_id].seq

        id_full = f'>{seq_id}.{Start}-{End}-{strand}-{gene_id}'
        # id_full_donor = f'>{seq_id}.{Start}-{End}-{strand}-{gene_id}-donor'
        # id_full_acceptor = f'>{seq_id}.{Start}-{End}-{strand}-{gene_id}-acceptor'

        if trimming_type == "splice2deep":
            with open(f'{filename_donor}_intron.fa', 'a') as donor_intron_file, open(f'{filename_acceptor}_intron.fa', 'a') as acceptor_intron_file:
                if args.flanking_size is not None:
                    flanking_size = args.flanking_size
                else:
                    flanking_size = flanking_size_dict[trimming_type]

                # start = Start - 300
                donor_start = Start - flanking_size
                acceptor_start = End - flanking_size + 2

                # end = End + 300
                donor_end = Start + flanking_size + 2
                acceptor_end = End + flanking_size

                if donor_start <= 0 or acceptor_start <= 0:
                    # start = 1
                    print(
                        f'Short start border on {id_full} {Start} {End} {contig_lengths[seq_id]}')
                else:
                    if acceptor_end >= contig_lengths[seq_id] or donor_end >= contig_lengths[seq_id]:
                        # end = contig_lengths[seq_id]
                        print(f'Short end border on {id_full}')
                    else:
                        # intron_seq = seq_full[start:end]
                        # donor_seq = seq_full[donor_start:donor_end]
                        # acceptor_seq = seq_full[acceptor_start:acceptor_end]

                        if strand == "+":
                            #intron_seq = intron_seq
                            donor_seq = seq_full[donor_start:donor_end].upper(
                            )
                            acceptor_start = acceptor_start - 4
                            acceptor_seq = seq_full[acceptor_start:acceptor_end].upper(
                            )
                        elif strand == "-":
                            acceptor_seq = seq_full[donor_start:donor_end].upper(
                            )
                            donor_seq = seq_full[acceptor_start:acceptor_end].upper(
                            )
                            donor_seq = Seq(donor_seq)
                            donor_seq = donor_seq.reverse_complement()
                            acceptor_seq = Seq(acceptor_seq)
                            acceptor_seq = acceptor_seq.reverse_complement()
                        else:
                            print(
                                f'Strand {strand} is not + or - (gene {gene_id} and\
                                        seq_id {seq_id}')

                        donor_intron_file.write(f'{id_full}\n')
                        donor_intron_file.write(f'{donor_seq}\n')

                        acceptor_intron_file.write(f'{id_full}\n')
                        acceptor_intron_file.write(f'{acceptor_seq}\n')
                        #print(f"\n{id_full} \n{donor_seq}\n{acceptor_seq}\n")
        elif trimming_type == "IMEter":
            # start = Start - 300

            if args.flanking_size is not None:
                flanking_size = args.flanking_size
            else:
                flanking_size = flanking_size_dict[trimming_type]

            intron_start = Start + flanking_size
            intron_end = End + flanking_size

            if intron_start <= 0:
                print(
                    f'Short start border on {id_full} {Start} {End} {contig_lengths[seq_id]}')
            else:
                if intron_end >= contig_lengths[seq_id]:
                    # end = contig_lengths[seq_id]
                    print(f'Short end border on {id_full}')
                else:
                    # intron_seq = seq_full[start:end]
                    # donor_seq = seq_full[donor_start:donor_end]
                    # acceptor_seq = seq_full[acceptor_start:acceptor_end]

                    intron_seq = seq_full[intron_start:intron_end].upper()

                    if strand == "-":
                        intron_seq = Seq(intron_seq)
                        intron_seq = intron_seq.reverse_complement()
                    elif strand != '+':
                        print(
                            f'Strand {strand} is not + or - (gene {gene_id} and\
                                seq_id {seq_id}')

                    intron_file.write(f'{id_full}\n')
                    intron_file.write(f'{intron_seq}\n')
                    #print(f"\n{id_full} \n{donor_seq}\n{acceptor_seq}\n")

        elif trimming_type == "spliceator":

            if args.flanking_size is not None:
                flanking_size = args.flanking_size
            else:
                flanking_size = flanking_size_dict[trimming_type]

            start = Start - flanking_size
            end = End + flanking_size

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
