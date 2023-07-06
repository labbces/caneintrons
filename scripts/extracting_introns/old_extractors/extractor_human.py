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
infile = args.input

with open(fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq_full = record.seq
        #print(seq_full)
        k = 0
        with open(infile, "r") as infile:
            for line in infile:
                line = line.replace('\n', '')
                k += 1
                print(k)
                numb, chr, start, end, as_type, numb2, strand = line.split(
                    '\t')
                filename = args.output + "_" + as_type
                Start = int(start)
                End = int(end)

                id_full = f'>{chr}.{Start}-{End}-{strand}-{as_type}'

                donor_start = Start - 69
                acceptor_start = End - 69
                donor_end = Start + 71
                acceptor_end = End + 71

                if strand == "+":
                    donor_seq = seq_full[donor_start:donor_end].upper()
                    acceptor_seq = seq_full[acceptor_start:acceptor_end].upper(
                    )
                else:
                    acceptor_seq = seq_full[donor_start:donor_end].upper()
                    donor_seq = seq_full[acceptor_start:acceptor_end].upper()
                    donor_seq = Seq(donor_seq)
                    donor_seq = donor_seq.reverse_complement()
                    acceptor_seq = Seq(acceptor_seq)
                    acceptor_seq = acceptor_seq.reverse_complement()
                    
                    #print(donor_seq)
                    #print(acceptor_seq)

                with open(f'{filename}_Donor_intron.fa', 'a') as donor_intron_file,  open(f'{filename}_Acceptor_intron.fa', 'a') as acceptor_intron_file:

                    donor_intron_file.write(f'{id_full}\n')
                    donor_intron_file.write(f'{donor_seq}\n')

                    acceptor_intron_file.write(f'{id_full}\n')
                    acceptor_intron_file.write(f'{acceptor_seq}\n')
