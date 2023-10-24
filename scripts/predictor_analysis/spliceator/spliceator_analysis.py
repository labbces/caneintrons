# Initial data: Spliceator output
# Objectives: 1) Getting the number of introns with 1 predicted splice site ok
#             2) Getting the number of introns with both (2) splice sites
#                predicted ("Good introns") OK
#             3) Getting the amount of genes with at least 1 "Good intron"
#             4) Getting the distribution of good introns per gene
# Output: Report of the four objectives

# imports
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import glob
import argparse
from Bio import motifs
from Bio.Seq import Seq
import pandas as pd
import logomaker

# Making it easier to load information through terminal
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--spliceator_file", type=str, required=True,
                    help='(Path to) Output from spliceator prediction')
parser.add_argument("-c", "--contig_table", type=str, required=True,
                    help='(Path to) table with contig length information in \
                        csv format. Ex: JXQF01000001.1,26708')
parser.add_argument("-e", "--extra", type=int, required=True,
                    help='Extra amount of nucleotides upstream/downstream \
                        intron donor splice site (e.g. 200, 50, ..)', default=200)
parser.add_argument("-d", "--plots_filename", type=str, required=True,
                    help='Filename of intron distribuition per gene histogram')
args = parser.parse_args()

# Loading initial data and defining variables
spliceator_f = args.spliceator_file

extra = args.extra
extra_up = extra - 5  # For 200 extra nucleotides, extra_up is 195
extra_down = extra + 5  # For 200 extra nucleotides, extra_down is 205

good_donors = {}
good_acceptors = {}
good_introns = {}
sense_sting_neg = {"don": [], "acce": []}
sense_sting_pos = {"don": [], "acce": []}
resume = {"pos": {"don": 0, "acce": 0}, "neg": {"don": 0, "acce": 0}}

contig_table = args.contig_table
contig_table_pd = pd.read_csv(contig_table, header=None)
contig_dict = dict(zip(contig_table_pd[0], contig_table_pd[1]))

distance = {'Donor': {}, 'Acceptor': {}}


def good_intron_info(seq, ss_type, intron_ID, pos_esperada, pos, strand, score):
    if intron_ID not in good_introns.keys():
        good_introns[intron_ID] = {'Donor': '', 'Acceptor': ''}

    pos_esperada = int(pos_esperada)
    pos = int(pos)

    if intron_ID not in distance[ss_type].keys():
        good_introns[intron_ID][ss_type] = [seq, pos, strand, score]
        distance[ss_type][intron_ID] = abs(pos_esperada - pos)
        # print(f'{intron_ID} {abs(pos_esperada - pos)} < {distance[ss_type][intron_ID]}')
    else:
        if abs(pos_esperada - pos) < distance[ss_type][intron_ID]:
            # print(f'{intron_ID} = {distance[ss_type][intron_ID]}')
            good_introns[intron_ID][ss_type] = [seq, pos, strand, score]
            distance[ss_type][intron_ID] = abs(pos_esperada - pos)
            # print(f'{intron_ID} = {distance[ss_type][intron_ID]}')

    if ss_type == "Donor":
        ss_types == 'Acceptor'
        if good_introns[intron_ID][ss_types] == '':
            good_introns[intron_ID][ss_types] = good_acceptors[intron_ID]

    elif ss_type == 'Acceptor':
        ss_types = 'Donor'
        if good_introns[intron_ID][ss_types] == '':
            good_introns[intron_ID][ss_types] = good_donors[intron_ID]


# Answering 4 objectives
# for filename in glob.glob('/home/bia/sugarcane_introns_local/1. analsing spliceator/*.spliceator.default.out'): # noqa: E501
with open(spliceator_f, 'r') as spliceator:
    # next(spliceator)  # skip header (e.g. first line)
    for line in spliceator:
        if not line.startswith('ID;#'):
            # Getting introns features from spliceator file
            fields = line.split(';')
            try:
                SS_type = fields[2]
            except:  # noqa: E722 E261
                continue
            seq = fields[3]
            position = int(fields[4])
            score = fields[5]
            score = score.replace('\n', '')
            score = float(score)
            ID = fields[0]

# Getting introns features from spliceator ID report
            m = re.match(
                r'([A-Za-z_.a0-9._]*)\.([0-9]*)-([0-9]*)-([\+\-])', ID)

            contig = m.group(1)
            intron_start = int(m.group(2))
            intron_end = int(m.group(3))
            strand = str((m.group(4)))

# Analysis for + strand
            if strand == "+":
                # Defining the correct amount of upstream and downstream nucleotides
                # to assess whether the splice site position is correct
                if intron_start >= extra:
                    upstream = extra
                else:
                    upstream = intron_start

                if intron_end + extra < contig_dict[contig]:
                    downstream = extra
                else:
                    downstream = contig_dict[contig] - intron_end

                # sequence total size: intron + extra upstream + extra downstream
                intron_containing_seqLength = upstream + \
                    ((intron_end - intron_start) + 1) + downstream

                if SS_type == "Donor":
                    resume['pos']["don"] += 1
                    sense_sting_pos['don'].append(int(position))
                    if intron_start > extra:  # non-border introns
                        if extra_up < int(position) < extra_down:
                            good_donors[ID] = [seq, position, strand, score]
                            if ID in good_acceptors:
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Donor", intron_ID=ID, pos=position, pos_esperada=extra, score=score)
                                # good_introns[ID] = 1
                    else:  # <= border-introns
                        if (intron_start - 5) < int(position) < (intron_start + 5):
                            good_donors[ID] = [seq, position, strand, score]
                            # good_intron[gene][ID]
                            # qtde good_intron len(good_intron)
                            # qtde de intron/gene len(gene)
                            # good_donor = [gene][ID] = 1
                            # len(good_intron[gene])(len good_donor[gene] + len good_acceptor[gene]) = genes com bas introns
                            if ID in good_acceptors:
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Donor", intron_ID=ID, pos=position, pos_esperada=intron_start, score=score)

                                # good_introns[ID] = 1

                elif SS_type == "Acceptor":
                    resume['pos']["acce"] += 1
                    sense_sting_pos['acce'].append(int(position))
                    if (contig_dict[contig] - intron_end) > extra:  # non-border introns
                        if (intron_containing_seqLength - extra_up) > int(position) > (intron_containing_seqLength - extra_down):
                            good_acceptors[ID] = [seq, position, strand, score]
                            if ID in good_donors:
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Acceptor", intron_ID=ID, pos=position, pos_esperada=extra, score=score)

                                # good_introns[ID] = 1
                    else:  # border introns
                        if (contig_dict[contig] - (intron_end - 5)) > int(position) > (contig_dict[contig] - (intron_end + 5)):
                            good_acceptors[ID] = [seq, position, strand, score]
                            if ID in good_donors:
                                k = (contig_dict[contig]) - intron_end
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Acceptor", intron_ID=ID, pos=position, pos_esperada=k, score=score)

                                # good_introns[ID] = 1

            if strand == "-":

                # print(strand)

                # Defining the correct amount of upstream and downstream nucleotides
                # to assess whether the splice site position is correct
                if intron_start >= extra:
                    upstream = extra
                else:
                    upstream = intron_start

                if intron_end + extra < contig_dict[contig]:
                    downstream = extra
                else:
                    downstream = contig_dict[contig] - intron_end

                # verificar dps que corrigir o extractor.py
                # sequence total size: intron + extra upstream + extra downstream
                intron_containing_seqLength = upstream + \
                    ((intron_end - intron_start) + 1) + downstream

                if SS_type == "Donor":
                    resume['neg']["don"] += 1
                    # print(f'Donor: {position}')
                    sense_sting_neg['don'].append(int(position))
                    if (intron_end + extra) < contig_dict[contig]:
                        if (extra_up) < int(position) < (extra_down):
                            good_donors[ID] = [seq, position, strand, score]
                            if ID in good_acceptors:
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Donor", intron_ID=ID, pos=position, pos_esperada=extra, score=score)

                                # good_introns[ID] = 1
                    else:  # border intron
                        if (contig_dict[contig] - (intron_end)) < int(position) < (contig_dict[contig] - (intron_end)):
                            good_donors[ID] = [seq, position, strand, score]
                            if ID in good_acceptors:
                                k = (contig_dict[contig]) - intron_end
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Donor", intron_ID=ID, pos=position, pos_esperada=k, score=score)

                if SS_type == "Acceptor":
                    resume['neg']["acce"] += 1
                    sense_sting_neg['acce'].append(int(position))
                    if intron_start > extra:  # non-border intron
                        if ((intron_containing_seqLength - extra_up)) >= int(position) >= ((intron_containing_seqLength - extra_down)):
                            good_acceptors[ID] = [seq, position, strand, score]
                            if ID in good_donors:
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Acceptor", intron_ID=ID, pos=position, pos_esperada=extra, score=score)

                    else:  # border intron
                        if (intron_start) < int(position) < (intron_start):
                            good_acceptors[ID] = [seq, position, strand, score]
                            if ID in good_donors:
                                # print('I have a pair!')
                                good_intron_info(
                                    seq=seq, strand=strand, ss_type="Acceptor", intron_ID=ID, pos=position, pos_esperada=intron_start, score=score)

                                # good_introns[ID] = 1
# Displaying results
# getting good introns amount
introns_count = len(good_introns.values())
print(introns_count)

# getting bad introns amount
bad_introns = (len(good_acceptors) + len(good_donors)) - introns_count
print(bad_introns)
# First displaying
# print(f'1. BAD INTRONS) Introns with one predicted site: {bad_introns}')
# print(f'2. GOOD INTRONS) Introns with both predicted sites: {introns_count}')
# print(f'3. GOOD INTRONS GENES) Amount of genes with at least 1good intron: {genes_good_introns}')
# print(f'3.1 Average of good introns per gene {introns_count/genes_good_introns}')

final_file = args.plots_filename + '_spliceator_analysis.csv'
introns_headers = []
with open(final_file, 'w') as file_final:
    file_final.write('ID,donor_sequence,acceptor_sequence,strand,mean_score')
    for intron_ID in good_introns.keys():
        full_header = (f'>{intron_ID}')
        if full_header not in introns_headers:
            sequence_a = good_introns[intron_ID]['Acceptor'][0]
            sequence_d = good_introns[intron_ID]['Donor'][0]
            strand = good_introns[intron_ID]['Acceptor'][2]
            score_a = good_introns[intron_ID]['Acceptor'][3]
            score_d = good_introns[intron_ID]['Donor'][3]
            mean_score = (score_a+score_d)/2

            file_final.write(
                f'\n{full_header},{sequence_d},{sequence_a},{strand},{mean_score}')

            if full_header in introns_headers:
                print(full_header)
            else:
                introns_headers.append(full_header)
        else:
            print('Check here')


# Predicted positions plot]
filename = args.plots_filename + 'PredictedPositions.png'
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

# Minus sense acceptor
acce_positions = sense_sting_neg['acce']
acce_positions_sorted = sorted(acce_positions)
axs[1, 1].hist(acce_positions_sorted, bins=10,
               range=(150, 400), color='lightpink')
axs[1, 1].set_xlabel('Posição do Sítio aceitor')
axs[1, 1].set_ylabel('Quantidade de sítios aceitores')
axs[1, 1].set_title('Aceitadores Fita - (~300-350)')
axs[1, 1].set_xticks(np.arange(150, 400, 25))
# plt.show()

# Minus sense donor
don_positions = sense_sting_neg['don']
don_positions_sorted = sorted(don_positions)
axs[1, 0].hist(don_positions_sorted, bins=10,
               range=(150, 400), color='lightblue')
axs[1, 0].set_xticks(range(150, 401, 25))
axs[1, 0].set_xlabel('Posição do Sítio doador')
axs[1, 0].set_ylabel('Quantidade de sítios doadores')
axs[1, 0].set_title('Doadores Fita - (~200)')
# .show()

# Plus sense acceptor
pos_positions = sense_sting_pos['acce']
pos_positions_sorted = sorted(pos_positions)
axs[0, 1].hist(pos_positions_sorted, bins=10,
               range=(150, 400), color='lightpink')
axs[0, 1].set_xticks(range(150, 401, 25))
axs[0, 1].set_xlabel('Posição do Sítio doador')
axs[0, 1].set_ylabel('Quantidade de sítios aceitadores')
axs[0, 1].set_title('Aceitadores Fita + (~ 300-350)')
# .show()

# PLus sense donnor
pos_positions = sense_sting_pos['don']
pos_positions_sorted2 = sorted(pos_positions)
axs[0, 0].hist(pos_positions_sorted2, bins=10,
               range=(150, 400), color='lightblue')
axs[0, 0].set_xticks(range(150, 401, 25))
axs[0, 0].set_xlabel('Posição do Sítio doador')
axs[0, 0].set_ylabel('Quantidade de sítios doadores')
axs[0, 0].set_title('Doadores Fita + (~200)')
# .show()

plt.tight_layout()
# plt.savefig(filename)
# plt.show()

plt.figure()
plt.hist([don_positions_sorted, acce_positions_sorted,
         pos_positions_sorted2, pos_positions_sorted],
         color=['lightblue', 'royalblue', 'lightsalmon', 'tomato'], bins=5, range=(150, 400),
         label=['Donors +', 'Acceptors +', 'Donors -', 'Acceptors -'])
plt.xticks(range(150, 401, 25))
plt.xlabel('Posições')
plt.ylabel('Quantidade de sítios doadores e aceitadores')
plt.legend()
plt.show()


donors = don_positions_sorted + pos_positions_sorted2
donors = sorted(donors)
acceptors = acce_positions_sorted + pos_positions_sorted
acceptors = sorted(acceptors)
plt.figure(figsize=(15, 10))
plt.hist([donors, acceptors],
         color=['royalblue', 'tomato'], bins=5, range=(150, 400),
         label=['Doadores', 'Aceitadores'])
plt.xticks(range(0, 551, 25))
plt.xlabel('Posições')
plt.ylabel('Quantidade de sítios doadores e aceitadores')
plt.legend()
plt.show()


# python3 -u Spliceator_analysis.py -f "/home/bia/sugarcane_introns_local/spliceator_anaysis/plus/*/*"  -n "/home/bia/sugarcane_introns_local/spliceator_anaysis/minus/*/*" -c /home/bia/sugarcane_introns_local/data/Contig_and_Strand_tables/CTBE_IQUSP_contig_table.csv -e 200 -d CTBE_IQUSP_GoodIntronsPerGene.png
