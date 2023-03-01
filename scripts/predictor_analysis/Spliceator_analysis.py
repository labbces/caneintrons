# TO-DO add spliceator position analysis plot resume

# Initial data: Spliceator output
# Objectives: 1) Getting the number of introns with 1 predicted splice site ok
#             2) Getting the number of introns with both (2) splice sites
#                predicted ("Good introns")
#             3) Getting the amount of genes with at least 1 "Good intron"
# Output: Report of the three objectives
#         Histogram of good introns amount per gene
#         Plot of spliceator predicted site to Acceptor and Donor data and by strand

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
                    help='(Path to) Output from spliceator prediction. Ex: "~/spliceator_outputs/Athaliana*"')
# parser.add_argument("-t", "--ss_table", type=str, required=True, help='(Path to)\
#  table with strand information in tsv format. Ex: -\tevm.TU.SCSP803280_000095463.2')
parser.add_argument("-c", "--contig_table", type=str, required=True,
                    help='(Path to) table with contig length information in \
                        csv format. Ex: JXQF01000001.1,26708')
parser.add_argument("-e", "--extra", type=int, required=True,
                    help='Extra amount of nucleotides upstream/downstream \
                        intron donor splice site (e.g. 200, 50, ..)')
parser.add_argument("-d", "--plots_filename", type=str, required=True,
                    help='Filename of intron distribuition per gene histogram. Ex: plots_Athaliana')
args = parser.parse_args()

# Loading files and defining variables
spliceator_output_files = args.spliceator_file

extra = args.extra
extra_up = extra - 5  # For 200 extra nucleotides, extra_up is 195
extra_down = extra + 5  # For 200 extra nucleotides, extra_down is 205

good_donors = {}   # Directory with all good donors intron ID
good_acceptors = {}  # Directory with all good acceptors intron ID
good_introns = {}  # Directory with all good introns gene IDs and intron ID

# Positions of donors (don) and acceptors (acce) in plus and minus strand and its amount
sense_sting_neg = {"don": [], "acce": []}
sense_sting_pos = {"don": [], "acce": []}
resume = {"pos": {"don": 0, "acce": 0}, "neg": {"don": 0, "acce": 0}}

contig_table = args.contig_table
contig_table_pd = pd.read_csv(contig_table, header=None)
contig_dict = dict(zip(contig_table_pd[0], contig_table_pd[1]))


# ANSWERING 3 OBJECTIVES
# Loading spliceator output files (Remind: It should be informed between "")
for filename in glob.glob(spliceator_output_files):
    with open(os.path.join(os.getcwd(), filename), 'r') as spliceator:
        for line in spliceator:
            if not line.startswith('ID;#'):  # skip header
                '''Getting introns features from spliceator file'''
                fields = line.split(';')  # all features
                SS_type = fields[2]  # Donor or acceptor
                seq = fields[3]
                position = fields[4]
                ID = fields[0]

                m = re.match(
                    r'([A-Za-z_.a0-9._]*)\.([0-9]*)-([0-9]*)-([\+\-])-(.*)', ID)

                contig = m.group(1)
                intron_start = int(m.group(2))
                intron_end = int(m.group(3))
                strand = str((m.group(4)))
                gene_ID = m.group(5)
                # gene_ID = m.group(5)
                #strand = ss_dict[gene_ID]

                '''Defining the correct amount of upstream and downstream nucleotides
                to assess whether the splice site position is correct'''

                # Identifying amount of nucleotides on borders
                if intron_start >= extra:
                    upstream = extra
                else:
                    upstream = intron_start

                if intron_end + extra < contig_dict[contig]:
                    downstream = extra
                else:
                    downstream = contig_dict[contig] - intron_end

                # Getting sequence total size: intron + extra upstream + extra downstream
                intron_containing_seqLength = upstream + \
                    ((intron_end - intron_start) + 1) + downstream

# Basic workfow:
# 1) separate by sequence strend (to analyze if the splice site
#    is in the correct position)
# 2) separate by splice site type (e.g. donor or acceptor) to
#    get intron with both predicted
# 3) check if splice site position is correct
#    3.1) Filling informational variables
#    3.2) for non-border introns (start coord > extra)
#    3.3) for border introns
# 4) correct sites are added to good donor/acceptors
# 5) if a good donor/acceptor is in good acceptor/donnor it is added to good introns

                '''**** PLUS STRAND ****'''
                if strand == "+":
                    if SS_type == "Donor":
                        # Filling info variables
                        resume['pos']["don"] += 1
                        sense_sting_pos['don'].append(int(position))

                        # Analyzing introns splice sites
                        if intron_start > extra:  # non-border introns
                            if extra_up < int(position) < extra_down:
                                good_donors[ID] = 1
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                        else:  # <= border-introns
                            if (intron_start - 5) < int(position) < (intron_start + 5):
                                good_donors[ID] = 1
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1

                    elif SS_type == "Acceptor":
                        # Filling info variables
                        resume['pos']["acce"] += 1
                        sense_sting_pos['acce'].append(int(position))

                        # Analyzing introns splice sites
                        if (contig_dict[contig] - intron_end) > extra:  # non-border introns
                            if (intron_containing_seqLength - extra_up) > int(position) > (intron_containing_seqLength - extra_down):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                        else:  # border introns
                            if (contig_dict[contig] - (intron_end - 5)) > int(position) > (contig_dict[contig] - (intron_end + 5)):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1

                '''**** MINUS STRAND****'''
                if strand == "-":
                    if SS_type == "Donor":
                        # Filling info variables
                        resume['neg']["don"] += 1
                        sense_sting_neg['don'].append(int(position))

                        # Analyzing introns splice sites
                        # non-border introns
                        if (intron_end + extra) < contig_dict[contig]:
                            if (extra_up) < int(position) < (extra_down):
                                good_donors[ID] = 1
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                        else:  # border intron
                            if (contig_dict[contig] - (intron_end)) < int(position) < (contig_dict[contig] - (intron_end)):
                                good_donors[ID] = 1
                                if ID in good_acceptors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1

                    if SS_type == "Acceptor":
                        # Filling info variables
                        resume['neg']["acce"] += 1
                        sense_sting_neg['acce'].append(int(position))

                        # Analyzing introns splice sites
                        if intron_start > extra:  # non-border intron
                            if ((intron_containing_seqLength - extra_up)) >= int(position) >= ((intron_containing_seqLength - extra_down)):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1
                        else:  # border intron
                            if (intron_start) < int(position) < (intron_start):
                                good_acceptors[ID] = 1
                                if ID in good_donors:
                                    if gene_ID not in good_introns:
                                        good_introns[gene_ID] = {}
                                    good_introns[gene_ID][ID] = 1

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

# Objectives displaying
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
    #print(f'{gene} {introns}')
    intron_amount = len(introns)
    if intron_amount not in hist:
        hist[intron_amount] = 1
    else:
        hist[intron_amount] += 1
print(hist)

'''*** PLOTS ***'''

# Histogram plotting and saving
filename = args.plots_filename + '_GoodIntronsPerGene.png'
plt.bar(hist.keys(), hist.values(), color='lightgreen')
plt.xlabel('Introns amount')
plt.ylabel('Gene amount')
plt.title('Distribuition of introns per gene on sugarcane')
# plt.savefig(filename)

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
plt.ylabel('Quantidade de sítios doadores')
plt.legend()
plt.show()


donors = don_positions_sorted + pos_positions_sorted2
donors = sorted(donors)
acceptors = acce_positions_sorted + pos_positions_sorted
acceptors = sorted(acceptors)
plt.figure()
plt.hist([donors, acceptors],
         color=['royalblue', 'tomato'], bins=5, range=(150, 400),
         label=['Donors', 'Acceptors'])
plt.xticks(range(150, 401, 25))
plt.xlabel('Posições')
plt.ylabel('Quantidade de sítios doadores')
plt.legend()
plt.show()
