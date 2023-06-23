from Bio import SeqIO
from Bio.Seq import Seq
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True,
                    help='File with Chromosome, start/end coordinates and ID')
parser.add_argument("-f", "--path_to_fasta", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

fasta = args.path_to_fasta
inx = args.path_to_fasta+'.genome_inx'
fasta_index = SeqIO.index_db(inx, fasta, 'fasta')
infile = args.input
# filename = args.output

with open(infile, 'r') as infile:
    a = []
    k = 0
    for line in infile:
        line = line.replace('\n', '')
        numb, chr, start, end, as_type, numb2, strand = line.split('\t')
        filename = args.output + "_" + as_type

        with open(f'{filename}_Donor_intron.fa', 'a') as donor_intron_file,  open(f'{filename}_Acceptor_intron.fa', 'a') as acceptor_intron_file:

            Start = int(start)
            End = int(end)

            if chr in fasta_index:
                print('Chromosome exists')
                seq_full = fasta_index[chr].seq
                # print(seq_full)
                id_full = f'>{chr}.{Start}-{End}-{strand}-{as_type}'

                donor_start = Start - 69
                acceptor_start = End - 69
                donor_end = Start + 71
                acceptor_end = End + 71

                if strand == '+' or strand == '-':
                    if strand == "+":
                        donor_seq = seq_full[donor_start:donor_end].upper()
                        acceptor_seq = seq_full[acceptor_start:acceptor_end].upper(
                        )
                    elif strand == "-":
                        acceptor_seq = seq_full[donor_start:donor_end].upper()
                        donor_seq = seq_full[acceptor_start:acceptor_end].upper(
                        )
                        # print(donor_seq)
                        # donor_seq = Seq(donor_seq)
                        donor_seq = donor_seq.reverse_complement()
                        # acceptor_seq = Seq(acceptor_seq)
                        acceptor_seq = acceptor_seq.reverse_complement()

                    donor_intron_file.write(f'{id_full}\n')
                    donor_intron_file.write(f'{donor_seq}\n')

                    acceptor_intron_file.write(f'{id_full}\n')
                    acceptor_intron_file.write(f'{acceptor_seq}\n')
                else:
                    print(f'Seq not + or - {id_full}')
