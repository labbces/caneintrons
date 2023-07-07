### COORDINATES PATTERNS ###
# For constitutive splicing (CS) - constitutive exon: chr1:4706-5095:[+/-] e.g. chromosome, start, end, strand
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
                    required=True, help="without chr")
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
                id, chr, c1, c2, strand, event = line.split(
                    '\t')  # c1 = coordenate 1
                chr = chr.upper()
                chr = chr.replace('CHR', '')
                if chr == cromo.upper():
                    print(f'Sequence {id} being processed...')
                    if "EX" in event:
                        with open(f'{filename}_EXsk_acceptor.fa', 'a') as eka, open(f'{filename}_EXsk_donor.fa', 'a') as ekd:
                            # print(full_co, strand)
                            c1 = int(float(c1))
                            c2 = int(float(c2))

                            # Defining fasta sequences headers
                            id_full = f'>{id}-{c1}-{c2}-{strand}-chr{chr}-{event}'

                            # Getting coordinates to extract sequences
                            exon_c1_start = c1 - 69
                            exon_c1_end = c1 + 71
                            C1_end_seq = extract(
                                seq_full, exon_c1_start, exon_c1_end, strand, 2)

                            if strand == "+":
                                eka.write(f'{id_full}\n')
                                eka.write(f'{C1_end_seq}\n')
                            else:  # if strand -
                                ekd.write(f'{id_full}\n')
                                ekd.write(f'{C1_end_seq}\n')

                            exon_c2_start = c2 - 69
                            exon_c2_end = c2 + 71
                            C2_end_seq = extract(
                                seq_full, exon_c2_start, exon_c2_end, strand, 2)

                            if strand == "+":
                                ekd.write(f'{id_full}\n')
                                ekd.write(f'{C1_end_seq}\n')
                            else:  # if strand -
                                eka.write(f'{id_full}\n')
                                eka.write(f'{C1_end_seq}\n')

                    if "INT" in event:
                        with open(f'{filename}_INT_donor.fa', 'a') as intd, open(f'{filename}_INT_acceptor.fa', 'a') as inta:
                            c1 = int(float(c1))
                            c2 = int(float(c2))

                            # Defining fasta sequences headers
                            id_full = f'>{id}.{c1}-{c2}-{strand}-chr{chr}-{event}'

                            # Getting coordinates to extract sequences
                            intron_c1_start = c1 - 69
                            intron_c1_end = c1 + 71
                            intron_c2_start = c2 - 69
                            intron_c2_end = c2 + 71

                            start_seq = extract(
                                seq_full, intron_c1_start, intron_c1_end, strand, 0)
                            if strand == "+":
                                intd.write(f'{id_full}\n')
                                intd.write(f'{start_seq}\n')
                            else:  # if strand -
                                inta.write(f'{id_full}\n')
                                inta.write(f'{start_seq}\n')

                            end_seq = extract(
                                seq_full, intron_c2_start, intron_c2_end, strand, 0)
                            if strand == "-":
                                intd.write(f'{id_full}\n')
                                intd.write(f'{end_seq}\n')
                            else:  # if strand +
                                inta.write(f'{id_full}\n')
                                inta.write(f'{end_seq}\n')

                    if "ALTD" in event:
                        with open(f'{filename}_ALTD_donor.fa', 'a') as atd:
                            id_full = f'>{id}.{c1}-{c2}-{strand}-chr{chr}-{event}'
                            if strand == "+":
                                c2 = int(float(c2))
                                extra_start = c2 - 69
                                extra_end = c2 + 71

                                seq_atd = extract(
                                    seq_full, extra_start, extra_end, strand, 0)

                                atd.write(f'{id_full}\n')
                                atd.write(f'{seq_atd}\n')
                            else:  # if strand "-"
                                c1 = int(float(c1))
                                extra_start = c1 - 69
                                extra_end = c1 + 71

                                seq_atd = extract(
                                    seq_full, extra_start, extra_end, strand, 0)

                                atd.write(f'{id_full}\n')
                                atd.write(f'{seq_atd}\n')

                    if "ALTA" in event:
                        with open(f'{filename}_ALTA_acceptor.fa', 'a') as ata:
                            id_full = f'>{id}.{c1}-{c2}-{strand}-chr{chr}-{event}'
                            if strand == "-":
                                c2 = int(float(c2))
                                extra_start = c2 - 69
                                extra_end = c2 + 71

                                seq_atd = extract(
                                    seq_full, extra_start, extra_end, strand, 0)

                                ata.write(f'{id_full}\n')
                                ata.write(f'{seq_atd}\n')
                            else:  # if strand "+"
                                c1 = int(float(c1))
                                extra_start = c1 - 69
                                extra_end = c1 + 71

                                seq_atd = extract(
                                    seq_full, extra_start, extra_end, strand, 0)

                                ata.write(f'{id_full}\n')
                                ata.write(f'{seq_atd}\n')

                    if "CS" in event:  # constitutive splicing
                        c1 = int(float(c1))
                        c2 = int(float(c2))
                        with open(f'{filename}_CS_donor.fa', 'a') as donor_intron_file,  open(f'{filename}_CS_acceptor.fa', 'a') as acceptor_intron_file:
                            if chr.upper() == cromo:
                                id_full = f'>chr{chr}.{c1}-{c2}-{strand}'

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
                                    # c2 because here we're working with exons, not introns
                                    donor_intron_file.write(f'{c2_seq}\n')

                                    acceptor_intron_file.write(f'{id_full}\n')
                                    acceptor_intron_file.write(
                                        f'{c1_seq}\n')  # exon

                                else:
                                    donor_intron_file.write(f'{id_full}\n')
                                    donor_intron_file.write(
                                        f'{c1_seq}\n')  # exon

                                    acceptor_intron_file.write(f'{id_full}\n')
                                    acceptor_intron_file.write(
                                        f'{c2_seq}\n')  # exon
