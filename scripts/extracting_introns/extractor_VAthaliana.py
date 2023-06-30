import re
from Bio import SeqIO
from Bio.Seq import Seq
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True,
                    help='pastDB')
parser.add_argument("-f", "--path_to_fasta", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
parser.add_argument("-c", "--chromosome", type=str,
                    required=True, help="no padrão do AS file")
args = parser.parse_args()

pastDB = args.input
#pastDB = '/home/bia/sugarcane_introns_local/data/Athaliana_genome/EVENT_INFO-araTha10.tab'
fasta = args.path_to_fasta
filename = args.output


def extract(seq_full, start, end, strand, neg_effect):
    if strand == "+":
        seq = seq_full[start:end].upper()
    elif strand == "-":
        seq = seq_full[(start + neg_effect):(end + neg_effect)].upper()
        seq = Seq(seq)
        seq = seq.reverse_complement()
    return seq


with open(fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq_full = record.seq
        with open(pastDB, "r") as infile:
            for line in infile:
                chr = line.split(sep='\t')[2]
                chr = chr.split(sep=':')[0].upper()
                if chr == args.chromosome.upper():
                    strand = line.split(sep='\t')[6]
                    strand = strand.split(sep=':')[2]
                    id = line.split('\t')[1]
                    print(f'Sequence {id} being processed...')
                    full_co = line.split('\t')[4]
                    alternative_co = line.split('\t')[2]

                    if "EX" in id:
                        with open(f'{filename}_EXrec_donor.fa', 'a') as erd, open(f'{filename}_EXsk_acceptor.fa', 'a') as eka, open(f'{filename}_EXsk_donor.fa', 'a') as ekd, open(f'{filename}_EXrec_acceptor.fa', 'a') as era:
                            m = re.match(
                                r'(chr.*):([0-9]*),([0-9+]*)-([0-9+]*),([0-9]*)', full_co)
                            if m:
                                # print(full_co, strand)
                                C1_end = int(m.group(2))
                                sk_start = m.group(3).split('+')
                                sk_start = [int(i) for i in sk_start]
                                sk_end = m.group(4).split('+')
                                sk_end = [int(i) for i in sk_end]
                                C2_start = int(m.group(5))

                                # Defining fasta sequences headers
                                id_full = f'>{id}.-{strand}-{chr}'

                                # Getting coordinates to extract sequences
                                C1_end_start = C1_end - 69
                                C1_end_end = C1_end + 71
                                C1_end_seq = extract(
                                    seq_full, C1_end_start, C1_end_end, strand, 2)
                                # print(C1_end_seq)
                                if strand == "+":
                                    erd.write(f'{id_full}\n')
                                    erd.write(f'{C1_end_seq}\n')
                                else:  # if strand -
                                    era.write(f'{id_full}\n')
                                    era.write(f'{C1_end_seq}\n')

                                C2_start_start = C1_end - 69
                                C2_start_end = C1_end + 71
                                C2_start_seq = extract(
                                    seq_full, C2_start_start, C2_start_end, strand, -3)
                                # print(C2_start_seq)
                                if strand == "-":
                                    erd.write(f'{id_full}\n')
                                    erd.write(f'{C1_end_seq}\n')
                                else:  # if strand +
                                    era.write(f'{id_full}\n')
                                    era.write(f'{C1_end_seq}\n')

                                for sk_strt in sk_start:
                                    sk_start_start = sk_strt - 69
                                    sk_start_end = sk_strt + 71
                                    sk_start_seq = extract(
                                        seq_full, sk_start_start, sk_start_end, strand, -3)
                                    # print(sk_start_seq)
                                    if strand == "-":
                                        ekd.write(f'{id_full}\n')
                                        ekd.write(f'{C1_end_seq}\n')
                                    else:  # if strand +
                                        eka.write(f'{id_full}\n')
                                        eka.write(f'{C1_end_seq}\n')

                                for sk_ed in sk_end:
                                    sk_end_start = sk_ed - 69
                                    sk_end_end = sk_ed + 71
                                    sk_end_seq = extract(
                                        seq_full, sk_end_start, sk_end_end, strand, 2)
                                    if strand == "+":
                                        ekd.write(f'{id_full}\n')
                                        ekd.write(f'{C1_end_seq}\n')
                                    else:  # if strand -
                                        eka.write(f'{id_full}\n')
                                        eka.write(f'{C1_end_seq}\n')

                    if "INT" in id:
                        with open(f'{filename}_INTdonor.fa', 'a') as intd, open(f'{filename}_INTacceptor.fa', 'a') as inta:
                            m = re.match(
                                r'(chr.*).*-([0-9]*)=([0-9]*)-.*:.*', full_co)
                            if m:
                                #print(full_co, strand)
                                start = int(m.group(2))
                                end = int(m.group(3))

                                # Defining fasta sequences headers
                                id_full = f'>{id}.{start}-{end}-{strand}-{chr}'

                                # Getting coordinates to extract sequences
                                start_start = start - 69
                                start_end = start + 71
                                end_start = end - 69
                                end_end = end + 71

                                start_seq = extract(
                                    seq_full, start_start, start_end, strand, 0)
                                if strand == "+":
                                    intd.write(f'{id_full}\n')
                                    intd.write(f'{start_seq}\n')
                                else:  # if strand -
                                    inta.write(f'{id_full}\n')
                                    inta.write(f'{start_seq}\n')

                                end_seq = extract(
                                    seq_full, end_start, end_end, strand, 0)
                                if strand == "-":
                                    intd.write(f'{id_full}\n')
                                    intd.write(f'{end_seq}\n')
                                else:  # if strand +
                                    inta.write(f'{id_full}\n')
                                    inta.write(f'{end_seq}\n')

                    if "ALTD" in id:
                        with open(f'{filename}_ALTdonor.fa', 'a') as atd:
                            alternative_co = alternative_co.split(":")[1]
                            start, end = alternative_co.split('-')

                            if strand == "+":
                                atstart = int(end) - 69
                                atend = int(end) + 71
                            else:  # if strand -
                                atstart = int(start) - 69
                                atend = int(start) + 71

                            seq_atd = extract(
                                seq_full, atstart, atend, strand, 0)
                            id_full = f'>{id}.{atstart}-{atend}-{strand}-{chr}'

                            atd.write(f'{id_full}\n')
                            atd.write(f'{seq_atd}\n')

                    if "ALTA" in id:
                        with open(f'{filename}_ALTacceptor.fa', 'a') as ata:
                            alternative_co = alternative_co.split(":")[1]
                            start, end = alternative_co.split('-')

                            if strand == "+":
                                atstart = int(start) - 69
                                atend = int(start) + 71
                            else:  # if strand -
                                atstart = int(end) - 69
                                atend = int(end) + 71

                            seq_ata = extract(
                                seq_full, atstart, atend, strand, 0)
                            id_full = f'>{id}.{atstart}-{atend}-{strand}-{chr}'

                            ata.write(f'{id_full}\n')
                            ata.write(f'{seq_ata}\n')
