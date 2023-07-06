### COORDINATES PATTERNS ###
# For constitutive splicing (CS): chr1:4706-5095:[+/-] e.g. chromosome, start, end, strand
# For exons skipping (EX): chr1:4605,4706-5095,5174:[+/-] e.g chromosome:recognized E1end,skipped E2start-skipped E2end,recognized E3start:strand
# For alternative donor (ALTD) - positive stand: chr2:-14152343:+ e.g chromosome:-alternative site:strand
# For alternative donor (ALTD) - negative stand: chr2:14152343-:- e.g chromosome:alternative site-:strand
# For alternative acceptor (ALTA) - positive stand: chr2:14152343-:+ e.g chromosome:alternative site-:strand
# For alternative accepor (ALTA) - negative stand: chr2:-14152343:- e.g chromosome:-alternative site:strand
# For intron retention (IR) - chr1:101010-202020:[+/-] e.g chromosome:intron start-intron end:strand

# Imports
import re
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

# argarse to make it easier to use on terminal
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True,
                    help='tsv file with Alternative_splicing_type\tcoordinates. For info coordinates patterns see script header')
parser.add_argument("-f", "--path_to_fasta", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
parser.add_argument("-c", "--chromosome", type=str,
                    required=True, help="no padrão do AS file")
args = parser.parse_args()

# Saving args into variables
input = args.input
fasta = args.path_to_fasta
filename = args.output
cromo = args.chromosome.upper()

# Defining extracting function


def extract(seq_full, start, end, strand, neg_effect):
    startt = min([start, end])
    endd = max([start, end])
    if strand == "+":
        seq = seq_full[startt:endd].upper()
    elif strand == "-":
        seq = seq_full[(startt + neg_effect):(endd + neg_effect)].upper()
        seq = Seq(seq)
        seq = seq.reverse_complement()
    return seq


# Openning fasta to get sequence to extract interest sub-sequences from it
with open(fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq_full = record.seq
        # Openning input to get coordinates to extract sub-sequences
        with open(input, "r") as infile:
            for line in infile:
                id = line.split(sep='\t')[0]
                coord = line.split(sep='\t')[1]
                chr = coord.split(sep=':')[0].upper()
                if chr == args.chromosome.upper():
                    print(f'Sequence {id} being processed...')

                    # strand = line.split(sep='\t')[6]
                    # strand = strand.split(sep=':')[2]
                    # id = line.split('\t')[1]
                    # full_co = line.split('\t')[4]
                    # alternative_co = line.split('\t')[2]

                    if "EX" in id:
                        with open(f'{filename}_EXrec_donor.fa', 'a') as erd, open(f'{filename}_EXsk_acceptor.fa', 'a') as eka, open(f'{filename}_EXsk_donor.fa', 'a') as ekd, open(f'{filename}_EXrec_acceptor.fa', 'a') as era:
                            m = re.match(
                                r'(.*):([0-9]*),([0-9+]*)-([0-9+]*),([0-9]*):([+/-])', coord)
                            if m:
                                # print(full_co, strand)
                                C1_end = int(m.group(2))
                                sk_start = m.group(3).split('+')
                                sk_start = [int(i) for i in sk_start]
                                sk_end = m.group(4).split('+')
                                sk_end = [int(i) for i in sk_end]
                                C2_start = int(m.group(5))
                                strand = m.group(6)

                                # Defining fasta sequences headers
                                id_full = f'>{id}-{coord.split(":")[1]}-{strand}-{chr}'

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
                        with open(f'{filename}_INT_donor.fa', 'a') as intd, open(f'{filename}_INT_acceptor.fa', 'a') as inta:
                            m = re.match(
                                r'(.*):([0-9]*)-([0-9]*):([+/-])', coord)
                            if m:
                                start = int(m.group(2))
                                end = int(m.group(3))
                                strand = m.group(4)

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
                        with open(f'{filename}_ALTD_donor.fa', 'a') as atd:
                            m = re.match(
                                r'(.*):.*,([0-9]*)-([0-9]*),([0-9]*):([+/-])', coord)
                            if m:
                                # alternative_co = alternative_co.split(":")[1]
                                # start, end = alternative_co.split('-')
                                c1 = m.group(2)
                                c2 = m.group(3)
                                c3 = m.group(4)
                                strand = m.group(5)

                                lista_ordenada = sorted(
                                    [c1, c2, c3], reverse=True)
                                c = lista_ordenada[1]

                                c_start = int(c) - 69
                                c_end = int(c) + 71

                                seq_atd = extract(
                                    seq_full, c_start, c_end, strand, 0)
                                id_full = f'>{id}.{c_start}-{c_end}-{strand}-{chr}'

                                atd.write(f'{id_full}\n')
                                atd.write(f'{seq_atd}\n')

                    if "ALTA" in id:
                        with open(f'{filename}_ALTA_acceptor.fa', 'a') as ata:
                            m = re.match(
                                r'(.*):([0-9]*),([0-9]*)-([0-9]*),.*:([+/-])', coord)
                            if m:
                                # alternative_co = alternative_co.split(":")[1]
                                # start, end = alternative_co.split('-')
                                c1 = m.group(2)
                                c2 = m.group(3)
                                c3 = m.group(4)
                                lista_ordenada = sorted(
                                    [c1, c2, c3], reverse=True)
                                c = lista_ordenada[1]
                                strand = m.group(5)
                                c_start = int(c) - 69
                                c_end = int(c) + 71

                            seq_ata = extract(
                                seq_full, c_start, c_end, strand, 0)
                            id_full = f'>{id}.{c_start}-{c_end}-{strand}-{chr}'

                            ata.write(f'{id_full}\n')
                            ata.write(f'{seq_ata}\n')

                    if "CS" in id:  # constitutive splicing
                        m = re.match(r'(.*):([0-9]*)-([0-9]*):([+/-]*)', coord)
                        chr = m.group(1)
                        with open(f'{filename}_CS_donor.fa', 'a') as donor_intron_file,  open(f'{filename}_CS_acceptor.fa', 'a') as acceptor_intron_file:
                            if chr.upper() == cromo:
                                c1 = m.group(2)
                                c2 = m.group(3)
                                strand = m.group(4)
                                c1 = int(c1)
                                c2 = int(c2)

                                id_full = f'>{chr}.{c1}-{c2}-{strand}'

                                c1_start = c1 - 69
                                c2_start = c2 - 69
                                c1_end = c1 + 71
                                c2_end = c2 + 71

                                c1_seq = extract(
                                    seq_full=seq_full, start=c1_start, end=c1_end, strand=strand, neg_effect=0)
                                c2_seq = extract(
                                    seq_full=seq_full, start=c2_start, end=c2_end, strand=strand, neg_effect=0)

                                if strand == "+":
                                    donor_intron_file.write(f'{id_full}\n')
                                    donor_intron_file.write(f'{c2_seq}\n')

                                    acceptor_intron_file.write(f'{id_full}\n')
                                    acceptor_intron_file.write(f'{c1_seq}\n')

                                else:
                                    donor_intron_file.write(f'{id_full}\n')
                                    donor_intron_file.write(f'{c1_seq}\n')

                                    acceptor_intron_file.write(f'{id_full}\n')
                                    acceptor_intron_file.write(f'{c2_seq}\n')
