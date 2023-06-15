from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import pandas as pd
import re


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help='File with Chromosome, start/end coordinates and ID')
parser.add_argument("-f", "--path_to_fasta", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

fasta = args.path_to_fasta
inx = args.path_to_fasta+'.genome_inx'
fasta_index = SeqIO.index_db(inx, fasta, 'fasta')
infile = args.input
filename = args.output

df = pd.read_csv(infile, sep='\t')
df = df[df['REF_CO'].notna()]
df = df['REF_CO']
print(df)

coord = df.to_dict()
coord = list(coord.values())

with open(f'{filename}_Acceptor_as.fa', 'w') as donor_intron_file,  open(f'{filename}_Donor_as.fa', 'w') as acceptor_intron_file:
    for coordenada in coord:
        m = re.match(r'(chr[0-9]*):[0-9]*,([0-9]*)-([0-9]*),[0-9]*:([+/-])', coordenada)
        chr = m.group(1)
        chr = chr.replace('c', 'C')
        Start = int(m.group(2))
        End = int(m.group(3))
        strand = m.group(4)
        #print(f'Chr {chr}, Start: {Start}, End: {End}')

     
        seq_full = fasta_index[chr].seq
        id_full = f'>{chr}.{Start}-{End}-{strand}'

        donor_start = Start - 35
        acceptor_start = End - 35
        donor_end = Start + 35
        acceptor_end = End + 35

        chr_len = len(seq_full)

        if donor_start <= 0 or acceptor_start <= 0:
            # start = 1
            print(f'Short start border on {id_full} {Start} {End}')
        else: 
            if acceptor_end >= chr_len or donor_end >= chr_len:
                print(f'Short end border on {id_full}')
            else:
                if strand == '+' or strand == '-':
                    if strand == "+":
                        donor_seq = seq_full[donor_start:donor_end].upper()
                        acceptor_seq = seq_full[acceptor_start:acceptor_end].upper()
                    elif strand == "-":
                        acceptor_seq = seq_full[donor_start:donor_end].upper()
                        donor_seq = seq_full[acceptor_start:acceptor_end].upper()
                        donor_seq = Seq(donor_seq)
                        donor_seq = donor_seq.reverse_complement()
                        acceptor_seq = Seq(acceptor_seq)
                        acceptor_seq = acceptor_seq.reverse_complement()

                    donor_intron_file.write(f'{id_full}')
                    donor_intron_file.write(f'{donor_seq}\n')

                    acceptor_intron_file.write(f'{id_full}')
                    acceptor_intron_file.write(f'{acceptor_seq}\n')
                else:
                    print(f'Seq not + or - {id_full}')