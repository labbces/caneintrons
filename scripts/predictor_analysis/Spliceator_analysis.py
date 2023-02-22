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
import pandas as pd

# Making it easier to load information through terminal
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--spliceator_file", type=str, required=True,
                    help='(Path to) Output from spliceator prediction')
parser.add_argument("-t", "--ss_table", type=str, required=True,
                    help='(Path to) table with strand information in tsv\
                         format. Ex: -\tevm.TU.SCSP803280_000095463.2')
parser.add_argument("-c", "--contig_table", type=str, required=True,
                    help='(Path to) table with contig length information in \
                        csv format. Ex: JXQF01000001.1,26708')
parser.add_argument("-e", "--extra", type=int, required=True,
                    help='Extra amount of nucleotides upstream/downstream \
                        intron donor splice site (e.g. 200, 50, ..)')
parser.add_argument("-d", "--histogram_filename", type=str, required=True,
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

ss_table = args.ss_table
ss_table_pd = pd.read_csv(ss_table, sep='\t', header=None)
ss_dict = dict(zip(ss_table_pd[1], ss_table_pd[0]))

contig_table = args.contig_table
contig_table_pd = pd.read_csv(contig_table, header=None)
contig_dict = dict(zip(contig_table_pd[0], contig_table_pd[1]))


# Answering 4 objectives
# for filename in glob.glob('/home/bia/sugarcane_introns_local/1. analsing spliceator/*.spliceator.default.out'):
for filename in glob.glob(spliceator_f):
    with open(os.path.join(os.getcwd(), filename), 'r') as spliceator:
        # next(spliceator)  # skip header (e.g. first line)
        for line in spliceator:
            if not line.startswith('ID;#'):

                # Getting introns features from spliceator file
                fields = line.split(';')  # all features
                SS_type = fields[2]  # Donor or acceptor
                seq = fields[3]
                position = fields[4]
                ID = fields[0]

# Getting introns features from spliceator ID report
                # m = re.match(r'([A-Za0-9.]*)\.([0-9]*)-([0-9]*)\
                # -([\+\-])-(.*)', ID)
                m = re.match(
                    r'([A-Za0-9.]*)\.([0-9]*)-([0-9]*)-([\+\-])-(.*)', ID)

                contig = m.group(1)
                intron_start = int(m.group(2))
                intron_end = int(m.group(3))
                # strand_id = (m.group(4))
                # gene_ID = m.group(5)
                gene_ID = m.group(5)
                strand = ss_dict[gene_ID]
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

# filtering introns by splice sites position (getting intron in right
# positions)
# and getting good introns (objective 2)
#
# Basic workfow:
# 1) separate by sequence strend (to analyze if the splice site
#    is in the correct position)
# 2) separate by splice site type (e.g. donor or acceptor) to
#    get intron with both predicted
# 3) check if splice site position is correct
#    3.1) for non-border introns (start coord > extra)
#    3.2) for border introns
# 4) correct sites are added to good donor/acceptors
# 5) if a good donor/acceptor is in good acceptor/donnor it is added to good introns

                if strand == "+":
                    if SS_type == "Donor":
                        resume['pos']["don"] += 1
                        sense_sting_pos['don'].append(int(position))
                        if intron_start > extra:  # non-border introns
                            if extra_up < int(position) < extra_down:
                                good_donors[ID] = 1
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1
                        else:  # <= border-introns
                            if (intron_start - 5) < int(position) < (intron_start + 5):
                                good_donors[ID] = 1
                                # good_intron[gene][ID]
                                # qtde good_intron len(good_intron)
                                # qtde de intron/gene len(gene)
                                # good_donor = [gene][ID] = 1
                                # len(good_intron[gene])(len good_donor[gene] + len good_acceptor[gene]) = genes com bas introns
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1

                    elif SS_type == "Acceptor":
                        resume['pos']["acce"] += 1
                        sense_sting_pos['acce'].append(int(position))
                        if (contig_dict[contig] - intron_end) > extra:  # non-border introns
                            if (intron_containing_seqLength - extra_up) > int(position) > (intron_containing_seqLength - extra_down):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1
                        else:  # border introns
                            if (contig_dict[contig] - (intron_end - 5)) > int(position) > (contig_dict[contig] - (intron_end + 5)):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1

                if strand == "-":
                    if SS_type == "Donor":
                        resume['neg']["don"] += 1
                        #print(f'Donor: {position}')
                        sense_sting_neg['don'].append(int(position))
                        if (intron_end + extra) < contig_dict[contig]:
                            if (extra_up - 10) < int(position) < (extra_down + 10):
                                good_donors[ID] = 1
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        #print('I have a pair!')
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1
                        else:  # border intron
                            if (contig_dict[contig] - (intron_end) + 10) < int(position) < (contig_dict[contig] - (intron_end) - 10):
                                good_donors[ID] = 1
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        #print('I have a pair!')
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1
                    if SS_type == "Acceptor":
                        resume['neg']["acce"] += 1
                        sense_sting_neg['acce'].append(int(position))
                        if intron_start > extra:  # non-border intron
                            if ((intron_containing_seqLength - extra_up) + 10) >= int(position) >= ((intron_containing_seqLength - extra_down) - 10):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    if gene_ID not in good_introns:
                                        #print('I have a pair!')
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1
                        else:  # border intron
                            if (intron_start - 10) < int(position) < (intron_start + 10):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    #print('I have a pair!')
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                                    # good_introns[ID] = 1

# print(good_acceptors)
# print(good_donors)
# print(len(good_introns))
# print(sense_sting_neg)
print(resume)

# Displaying results
# getting good introns amount
introns_count = 0
for gene in good_introns.values():
    introns_count += len(gene)

# getting bad introns amount
bad_introns = (len(good_acceptors) + len(good_donors)) - introns_count


# Getting amount of genes with good introns
genes_good_introns = len(good_introns)

# First displaying
print(f'1. BAD INTRONS) Introns with one predicted site: {bad_introns}')
print(f'2. GOOD INTRONS) Introns with both predicted sites: {introns_count}')
print(
    f'3. GOOD INTRONS GENES) Amount of genes with at least 1\
    good intron: {genes_good_introns}')
print(
    f'3.1 Average of good introns per gene {introns_count/genes_good_introns}')


# Getting distribution hystogram
hist = {}
for gene, introns in good_introns.items():
    # print(f'{gene} {introns}')
    intron_amount = len(introns)
    if intron_amount not in hist:
        hist[intron_amount] = 1
    else:
        hist[intron_amount] += 1
# print(hist)

# Histogram plotting and saving
filename = args.histogram_filename
plt.bar(hist.keys(), hist.values(), color='lightgreen')
plt.xlabel('Introns amount')
plt.ylabel('Gene amount')
plt.title('Distribuition of introns per gene on sugarcane')
plt.savefig(filename)


fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

acce_positions = sense_sting_neg['acce']
acce_positions_sorted = sorted(acce_positions)
axs[1, 1].hist(acce_positions_sorted, bins=10,
               range=(150, 400), color='lightpink')
axs[1, 1].set_xlabel('Posição do Sítio aceitor')
axs[1, 1].set_ylabel('Quantidade de sítios aceitores')
axs[1, 1].set_title('Aceitadores Fita - (~300-350)')
axs[1, 1].set_xticks(np.arange(150, 400, 25))
# plt.show()

don_positions = sense_sting_neg['don']
don_positions_sorted = sorted(don_positions)
axs[1, 0].hist(don_positions_sorted, bins=10,
               range=(150, 400), color='lightblue')
axs[1, 0].set_xticks(range(150, 401, 25))
axs[1, 0].set_xlabel('Posição do Sítio doador')
axs[1, 0].set_ylabel('Quantidade de sítios doadores')
axs[1, 0].set_title('Doadores Fita - (~200)')
# .show()

pos_positions = sense_sting_pos['acce']
pos_positions_sorted = sorted(pos_positions)
axs[0, 1].hist(pos_positions_sorted, bins=10,
               range=(150, 400), color='lightpink')
axs[0, 1].set_xticks(range(150, 401, 25))
axs[0, 1].set_xlabel('Posição do Sítio doador')
axs[0, 1].set_ylabel('Quantidade de sítios aceitadores')
axs[0, 1].set_title('Aceitadores Fita + (~ 300-350)')
# .show()

pos_positions = sense_sting_pos['don']
pos_positions_sorted = sorted(pos_positions)
axs[0, 0].hist(pos_positions_sorted, bins=10,
               range=(150, 400), color='lightblue')
axs[0, 0].set_xticks(range(150, 401, 25))
axs[0, 0].set_xlabel('Posição do Sítio doador')
axs[0, 0].set_ylabel('Quantidade de sítios doadores')
axs[0, 0].set_title('Doadores Fita + (~200)')
# .show()

plt.tight_layout()
plt.show()
