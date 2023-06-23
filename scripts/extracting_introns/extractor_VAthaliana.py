import re
from Bio import SeqIO
from Bio.Seq import Seq

fasta = '/home/bia/sugarcane_introns_local/data/Genomas/Athaliana_whole_genome.fa'
pastDB = '/home/bia/sugarcane_introns_local/data/Athaliana_genome/EVENT_INFO-araTha10.tab'

inx = fasta+'.genome_inx'
fasta_index = SeqIO.index_db(inx, fasta, 'fasta')


def extract(seq_full, start, end, strand, neg_effect):
    if strand == "+":
        seq = seq_full[start:end].upper()
    elif strand == "-":
        seq = seq_full[(start + neg_effect):(end + neg_effect)].upper()
        seq = Seq(seq)
        seq = seq.reverse_complement()
    return seq


with open(pastDB, 'r') as pastDB, open('outfile.fa', 'w') as outfasta:
    with open('INT_donor.fa', 'w') as intD, open('INT_acceptor.fa', 'w') as intA:
        for line in pastDB:
            id = line.split('\t')[1]
            full_co = line.split('\t')[4]
            '''if "EX" in id:
                m = re.match(r'(chr[0-9]*):([0-9]*),([0-9+]*)-([0-9+]*),([0-9]*)', full_co)
                if m:
                    chr = m.group(1)
                    chr = chr.replace('c', 'C')
                    rec_donor = int(m.group(2))
                    sk_acceptor = m.group(3).split('+')
                    sk_acceptor = [int(i) for i in sk_acceptor]
                    sk_donor = m.group(4).split('+')
                    sk_donor = [int(i) for i in sk_donor]
                    rec_acceptor = int(m.group(5))

                    # getting chromosome sequence
                    seq_full = fasta_index[chr].seq
                    # print(seq_full)
                    
                    # Defining sequence strand
                    if rec_donor < rec_acceptor:
                        strand = '+'
                    else:
                        strand = '-'

                    # Defining fasta sequences headers
                    id_full = f'>{id}.-{strand}-{chr}'

                    # Getting coordinates to extract sequences
                    rec_donor_start = rec_donor - 34
                    rec_donor_end = rec_donor + 36
                    rec_donor_seq = extract(seq_full, rec_donor_start, rec_donor_end, strand, 2)
                    #print(rec_donor_seq)

                    #intD.write(f'{id_full}\n')
                    #intD.write(f'{rec_donor_seq}\n')


                    rec_acceptor_start = rec_donor - 34
                    rec_acceptor_end = rec_donor + 36
                    rec_acceptor_seq = extract(seq_full, rec_acceptor_start, rec_acceptor_end, strand, -3)
                    #print(rec_acceptor_seq)

                    #intA.write(f'{id_full}\n')
                    #intA.write(f'{rec_acceptor_seq}\n')

                    for acceptor in sk_acceptor:
                        sk_acceptor_start = acceptor - 34
                        sk_acceptor_end = acceptor + 36
                        sk_acceptor_seq = extract(seq_full, rec_acceptor_start, rec_acceptor_end, strand, -3)

                        #intA.write(f'{id_full}\n')
                        #intA.write(f'{sk_acceptor_seq}\n')

                    for donor in sk_donor:
                        sk_donor_start = donor - 34
                        sk_donor_end = donor + 36
                        sk_donor_seq = extract(seq_full, sk_donor_start, sk_donor_end, strand, 2)

                        #print(sk_donor_seq)

                        #intD.write(f'{id_full}\n')
                        #intD.write(f'{sk_donor_seq}\n')


            if "INT" in id:
                m = re.match(r'(chr[0-9]*).*-([0-9]*)=([0-9]*)-.*:.*', full_co)
                if m:
                    chr = m.group(1)
                    chr = chr.replace('c', 'C')
                    start = int(m.group(2))
                    end = int(m.group(3))

                    # getting chromosome sequence
                    seq_full = fasta_index[chr].seq
                    # print(seq_full)
                    
                    # Defining sequence strand
                    if start < end:
                        strand = '+'
                    else:
                        strand = '-'

                    # Defining fasta sequences headers
                    id_full = f'>{id}.{start}-{end}-{strand}-{chr}'

                    # Getting coordinates to extract sequences
                    donor_start = start - 34
                    donor_end =  start + 36

                    acceptor_start = end - 34
                    acceptor_end =  end + 36
                    
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
                        # print(acceptor_seq)
                    
                    #intD.write(f'{id_full}\n')
                    #intD.write(f'{donor_seq}\n')

                    #intA.write(f'{id_full}\n')
                    #intA.write(f'{acceptor_seq}\n')'''

            if "ALTD" in id:
                m = re.match(
                    r'(chr[0-9]*):[0-9+]*-([0-9+]*),([0-9]*)', full_co)
                if m:
                    chr = m.group(1)
                    chr = chr.replace('c', 'C')
                    start = m.group(2).split("+")
                    start = [int(i) for i in start]
                    acceptor = int(m.group(3))

                    # getting chromosome sequence
                    seq_full = fasta_index[chr].seq
                    # print(seq_full)

                    # Defining sequence strand
                    if start[0] < acceptor:
                        strand = '+'
                    else:
                        strand = '-'

                    # Defining fasta sequences headers
                    id_full = f'>{id}.{start}-{strand}-{chr}_ALTD'

                    # Getting coordinates to extract sequences
                    donor_start = start - 34
                    donor_end = start + 36

                    if strand == "+":
                        donor_seq = seq_full[donor_start:donor_end].upper()
                    elif strand == "-":
                        acceptor_seq = seq_full[donor_start:donor_end].upper()
                        acceptor_seq = Seq(acceptor_seq)
                        acceptor_seq = acceptor_seq.reverse_complement()
                        print(acceptor_seq)

                    # intD.write(f'{id_full}\n')
                    # intD.write(f'{donor_seq}\n')

                    # intA.write(f'{id_full}\n')
                    # intA.write(f'{acceptor_seq}\n')
