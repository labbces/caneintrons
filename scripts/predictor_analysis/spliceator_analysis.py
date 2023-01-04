# Initial data: Spliceator output
# Objectives: ok 1) Getting the number of introns with 1 predicted splice site
#             ok 2) Getting the number of introns with both (2) splice sites
#                predicted ("Good introns")
#              ok 3) Getting the distribution of good intros per gene
#             ok 4) Getting the amount of genes with at least 1 "Good intron"
# Output: Report of the four objectives

# imports
import matplotlib.pyplot as plt
import re
import os
import glob
import argparse

# Making it easier to load information through terminal
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--spliceator_file", type=str, required=True,
                    help='(Path to) Output from spliceator prediction')
parser.add_argument("-e", "--upstream_nt_extra",
                    type=int, required=True, help='Extra amount of nucleotides upstream intron donor splice site (e.g. 200, 50, ..)')
parser.add_argument("-d", "--histogram_filename", type=str, required=True,
                    help='Filename of intron distribuition per gene histogram')
args = parser.parse_args()

# Loading initial data
spliceator_f = args.spliceator_file
data = ''
# for filename in glob.glob('/home/bia/sugarcane_introns_local/1. analsing spliceator/*.spliceator.default.out'):
for filename in glob.glob(spliceator_f):
    with open(os.path.join(os.getcwd(), filename), 'r') as spliceator:
        next(spliceator)  # skip header (e.g. first line)
        spliceator_lines = spliceator.readlines()
        if data != '':
            data.extend(spliceator_lines)
        else:
            data = list(spliceator_lines)

# Objectives 1 (and 2)

# Separating acceptors and donors lines
acceptors = [line for line in data if "Acceptor" in line]
donors = [line for line in data if "Donor" in line]

# Getting good and bad introns
good_introns = [line for line in data if line.split(
    ';')[0] in '\t'.join(acceptors) and '\t'.join(donors)]


def trim(line):
    line = line.split(';')[0]
    return line


bad_introns = [*set([trim(line) for line in data if line not in good_introns])]
print('1. Introns with 1 predicted site: ',
      len(bad_introns))  # 1 objective result

# Objective 2
# Analysing good introns
ids = []
distrib = {}
splice_sites = {'Acceptor': {}, 'Donor': {}}
extra = args.upstream_nt_extra
extra_down = extra - 3
extra_up = extra + 3
for line in good_introns:
    m = re.match(r'(.*);#.*;(.*);(.*);(.*);.*', line)
    id = m.group(1)
    gene = id.split('_')[1]
    SS_type = m.group(2)
    sequence = m.group(3)
    position = m.group(4)
    site = sequence[10:12]

    if site in splice_sites[SS_type].keys():
        splice_sites[SS_type][site] += 1
    else:
        splice_sites[SS_type][site] = 1
    # print(f'{gene} {id} {SS_type} {sequence} {position}')

    if SS_type == "Donor":
        if extra_up > int(position) >= extra_down:
            if id not in ids:
                ids.append(id)
        else:
            print(
                f'Intron {id} donor splice site in a wrong position {position}')

    if gene not in distrib:
        distrib[gene] = [id]
    else:
        if id not in distrib[gene]:
            distrib[gene].append(id)
# 2 objective result
print(f'2. Introns with both predicted sites: {len(ids)}')

# Objective 3 (and 4)

for gene in distrib:
    distrib[gene] = len(distrib[gene])
distrib2 = {}
for amount in distrib.values():
    if amount not in distrib2:
        distrib2[amount] = 1
    else:
        distrib2[amount] += 1
print(
    f'3. The distribuition is (intron amount: gene amount) -see histogram for visual representation-: {distrib2}')  # Objective 4 result

distrib2 = [key for key, val in distrib2.items() for _ in range(val)]

filename = args.histogram_filename
plt.hist(distrib2, color='lightpink', align='left')
plt.xlabel('Introns amount')
plt.ylabel('Gene amount')
plt.title('Distribuition of introns per gene on sugarcane')
plt.savefig(filename)  # Objective 3 histogram

# Objective 4

print(
    f'4. The total of genes with at least 1 good intron is: {len((distrib).keys())}')  # Objective 4 result
print(f'Extra: introns splice sites frequency: {splice_sites}')

