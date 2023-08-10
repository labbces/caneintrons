### COORDINATES PATTERNS ###
# all coordinates files should have 6 columns: ID, CHR, C1, C2, Strand, Event
# in this specific order and without header
# C1 and C2 comprehend coordinates 1 and 2, they do not need to be in 5' - 3' order
# Strand should be either + or -
# Event can be: INT for reteined intron, ALTA for alternative acceptor/5', ALTD for alternative donor/3', EX for exon skipping or CS for constitutive splicing
# All lines should have a value, except for ALTA/ALTD events:
# For ALTA strand +, C2 can be empty | for ALTA strand -, C1 can be empty
# For ALTD strand +, C1 can be empty | for ALTD strand -, C2 can be empty
# For CS events the script expects an EXON coordinate
# For INT events the script expects an INTRON coordinate
# For EX events the script expects an EXONS coordinate
# For ALTA and ALTD the scrips expect the alternative SS site
# For an illustrative representation see figure "splicing_events.png" at projects GitHub
# We recommend to use the same pattern for chr name in gtf, genome and coordinate files
# ! CHR column must have the same name as fasta headers


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
parser.add_argument("-t", "--type", default='fasta_iterator', choices=['fasta_iterator', 'fasta_index'],
                    help="Way of selecting sequences from fasta file. 'fasta_iterator' is a better option for  fasta files with long sequences, such as scaffold or chromosome levels. For fasta files with smaller sequences such as contig sequences, 'fasta_index' is the best option - it will take some time to generate the index, but after that it is usualy faster. Default is 'fasta_interator'")
args = parser.parse_args()

# Saving args into variables
input = args.input
fasta = args.path_to_fasta
filename = args.output
cromo = args.chromosome.upper()

tipo = args.type

if tipo == "fasta_index":
    inx = args.path_to_fasta+'.genome_inx'
    fasta_index = SeqIO.index_db(inx, fasta, 'fasta')
    print("index done")


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


def full_seq(tipo, chr, fasta):
    output = []
    if tipo == "fasta_iterator":
        with open(fasta) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                header = record.id
                if header == chr:
                    seq_full = record.seq
                    output.append(header)
                    output.append(seq_full)
                    return output
    elif tipo == "fasta_index":
        seq_full = fasta_index[chr].seq
        header = chr
        output.append(header)
        output.append(seq_full)
        return output


# Openning fasta to get sequence to extract interest sub-sequences from it
with open(input, "r") as infile:
    for line in infile:
        id, chr, c1, c2, strand, event = line.rstrip('\n').split(
            '\t')  # c1 = coordenate 1
        chr = chr.strip()

        header = full_seq(tipo=tipo, chr=chr, fasta=fasta)[0]
        seq_full = full_seq(tipo=tipo, chr=chr, fasta=fasta)[1]
        # print(header, seq_full)
        seq_size = len(seq_full)

        id = id.strip()
        c1 = c1.strip()
        c2 = c2.strip()
        strand = strand.strip()
        event = event.strip()
        print(f'Sequence {id} being processed...')

        if "EX" in event:
            with open(f'{filename}_EXsk_acceptor.fa', 'a') as eka, open(f'{filename}_EXsk_donor.fa', 'a') as ekd:
                # print(full_co, strand)
                c1 = int(float(c1))
                c2 = int(float(c2))

                # Defining fasta sequences headers
                id_full = f'>{id}.{event}_{chr}:{c1}-{c2}:{strand}'

                # Getting coordinates to extract sequences
                exon_c1_start = c1 - 69
                exon_c1_end = c1 + 71

                if exon_c1_start > 0 and exon_c1_end > 0:
                    if exon_c1_start < seq_size and exon_c1_end < seq_size:
                        C1_end_seq = extract(
                            seq_full, exon_c1_start, exon_c1_end, strand, 2)

                        if strand == "+":
                            eka.write(f'{id_full.strip()}\n')
                            eka.write(f'{C1_end_seq.strip()}\n')
                        else:  # if strand -
                            ekd.write(f'{id_full.strip()}\n')
                            ekd.write(f'{C1_end_seq.strip()}\n')

                exon_c2_start = c2 - 69
                exon_c2_end = c2 + 71
                if exon_c2_start > 0 and exon_c2_end > 0:
                    if exon_c2_start < seq_size and exon_c2_end < seq_size:
                        C2_end_seq = extract(
                            seq_full, exon_c2_start, exon_c2_end, strand, 2)

                        if strand == "+":
                            ekd.write(f'{id_full.strip()}\n')
                            ekd.write(f'{C2_end_seq.strip()}\n')
                        else:  # if strand -
                            eka.write(f'{id_full.strip()}\n')
                            eka.write(f'{C2_end_seq.strip()}\n')

        if "INT" in event or 'INTRON' in event:
            with open(f'{filename}_INT_donor.fa', 'a') as intd, open(f'{filename}_INT_acceptor.fa', 'a') as inta:
                c1 = int(float(c1))
                c2 = int(float(c2))

                # Defining fasta sequences headers
                id_full = f'>{id}.{event}_{chr}:{c1}-{c2}:{strand}'

                # Getting coordinates to extract sequences
                intron_c1_start = c1 - 69
                intron_c1_end = c1 + 71
                intron_c2_start = c2 - 69
                intron_c2_end = c2 + 71

                if intron_c1_start > 0 and intron_c1_end > 0:
                    if intron_c1_start < seq_size and intron_c1_end < seq_size:
                        start_seq = extract(
                            seq_full, intron_c1_start, intron_c1_end, strand, 0)
                        if strand == "+":
                            intd.write(f'{id_full.strip()}\n')
                            intd.write(f'{start_seq.strip()}\n')
                        else:  # if strand -
                            inta.write(f'{id_full.strip()}\n')
                            inta.write(f'{start_seq.strip()}\n')

                if intron_c2_start > 0 and intron_c2_start > 0:
                    if intron_c2_start < seq_size and intron_c2_start < seq_size:
                        end_seq = extract(
                            seq_full, intron_c2_start, intron_c2_end, strand, 0)
                        if strand == "-":
                            intd.write(f'{id_full.strip()}\n')
                            intd.write(f'{end_seq.strip()}\n')
                        else:  # if strand +
                            inta.write(f'{id_full.strip()}\n')
                            inta.write(f'{end_seq.strip()}\n')

        if "ALTD" in event:
            with open(f'{filename}_ALTD_donor.fa', 'a') as atd:
                id_full = f'>{id}.{event}_{chr}:{c1}-{c2}:{strand}'
                if strand == "+":
                    c2 = int(float(c2))
                    extra_start = c2 - 69
                    extra_end = c2 + 71

                    if extra_start > 0 and extra_end > 0:
                        if extra_start < seq_size and extra_end < seq_size:
                            seq_atd = extract(
                                seq_full, extra_start, extra_end, strand, 0)
                            atd.write(f'{id_full.strip()}\n')
                            atd.write(f'{seq_atd.strip()}\n')

                else:  # if strand "-"
                    c1 = int(float(c1))
                    extra_start = c1 - 69
                    extra_end = c1 + 71

                    if extra_start > 0 and extra_end > 0:
                        if extra_start < seq_size and extra_end < seq_size:
                            seq_atd = extract(
                                seq_full, extra_start, extra_end, strand, 0)

                            atd.write(f'{id_full.strip()}\n')
                            atd.write(f'{seq_atd.strip()}\n')

        if "ALTA" in event:
            with open(f'{filename}_ALTA_acceptor.fa', 'a') as ata:
                id_full = f'>{id}.{event}_{chr}:{c1}-{c2}:{strand}'
                if strand == "-":
                    c2 = int(float(c2))
                    extra_start = c2 - 69
                    extra_end = c2 + 71

                    if extra_start > 0 and extra_end > 0:
                        if extra_start < seq_size and extra_end < seq_size:
                            seq_atd = extract(
                                seq_full, extra_start, extra_end, strand, 0)

                            ata.write(f'{id_full.strip()}\n')
                            ata.write(f'{seq_atd.strip()}\n')
                else:  # if strand "+"
                    c1 = int(float(c1))
                    extra_start = c1 - 69
                    extra_end = c1 + 71
                    if extra_start > 0 and extra_end > 0:
                        if extra_start < seq_size and extra_end < seq_size:
                            seq_atd = extract(
                                seq_full, extra_start, extra_end, strand, 0)

                            ata.write(f'{id_full.strip()}\n')
                            ata.write(f'{seq_atd.strip()}\n')

        if "CS" in event or "EXON in event":  # constitutive splicing
            c1 = int(float(c1))
            c2 = int(float(c2))
            with open(f'{filename}_CS_donor.fa', 'a') as donor_intron_file,  open(f'{filename}_CS_acceptor.fa', 'a') as acceptor_intron_file:
                if chr.upper() == cromo:
                    id_full = f'>{id}.{event}_{chr}:{c1}-{c2}:{strand}'

                    c1_start = c1 - 69
                    c2_start = c2 - 69
                    c1_end = c1 + 71
                    c2_end = c2 + 71

                    c1_seq = extract(
                        seq_full=seq_full, start=c1_start, end=c1_end, strand=strand, neg_effect=0)
                    c2_seq = extract(
                        seq_full=seq_full, start=c2_start, end=c2_end, strand=strand, neg_effect=0)

                    if strand == "+":
                        if c2_start > 0 and c2_end > 0:
                            if c2_start < seq_size and c2_end < seq_size:
                                donor_intron_file.write(
                                    f'{id_full.strip()}\n')
                                # c2 because here we're working with exons, not introns
                                donor_intron_file.write(
                                    f'{c2_seq.strip()}\n')

                        if c1_start > 0 and c1_end > 0:
                            if c1_start < seq_size and c1_end < seq_size:
                                acceptor_intron_file.write(
                                    f'{id_full.strip()}\n')
                                acceptor_intron_file.write(
                                    f'{c1_seq.strip()}\n')  # exon

                    else:
                        if c1_start > 0 and c1_end > 0:
                            if c1_start < seq_size and c1_end < seq_size:
                                donor_intron_file.write(
                                    f'{id_full.strip()}\n')
                                donor_intron_file.write(
                                    f'{c1_seq.strip()}\n')  # exon

                        if c2_start > 0 and c2_end > 0:
                            if c2_start < seq_size and c2_end < seq_size:
                                acceptor_intron_file.write(
                                    f'{id_full.strip()}\n')
                                acceptor_intron_file.write(
                                    f'{c2_seq.strip()}\n')  # exon
