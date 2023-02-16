# Initial data: RNAseq snakemake output - SJ.out.tab file
# Objectives: Get expression of junctions and report sites with
#             at least 20 supporting RNASeq reads
# Basic Workflow:
#   1) Load all SJ files
#   2) Get expression
#   3) Final report: Contig-Sites-Strend-ReadsNumber
# Output: Report of the four objectives
#
# How to run: python3 /home/user/RNAseq_analysis.py -p '/home/user/*SJ.out.tab' -o output.txt
#


# imports
import argparse
import pandas as pd
import glob
import os

# Making it easier to load information through terminal
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path_to_SJfile", type=str, required=True,
                    help='Path to output of STAR SJ file')
parser.add_argument("-o", "--out_prefixName", type=str, required=True,
                    help='Prefix name to output')
args = parser.parse_args()

SJ_files = args.path_to_SJfile

# Creating important variables
reads = {}
good_reads = {}

# Getting expression info from tab file to python dict
for filename in glob.glob(SJ_files):
    with open(os.path.join(os.getcwd(), filename), 'r') as SJ_file:
        sj_table_pd = pd.read_csv(SJ_file, sep='\t', header=None)
        sj_table_dict = sj_table_pd.to_dict()

        # Getting expression and well supported splice junction (at least 20 RNAseq mapping)
        n = 0
        while n < len(sj_table_dict[0]):
            # Creating id to save and linking it to its RNAseq mapping
            sj_id = f'{sj_table_dict[0][n]}_{((sj_table_dict[1][n])-1)}_{((sj_table_dict[2][n])-1)}_{sj_table_dict[3][n]}'

            uniquely_mapped = sj_table_dict[6][n]
            multi_mapped = sj_table_dict[7][n]
            mapped_sum = uniquely_mapped + multi_mapped

            # getting good SpliceJunctions and all info
            if mapped_sum > 0:
                if mapped_sum < 20:
                    if sj_id not in good_reads:
                        if sj_id not in reads.keys():
                            reads[sj_id] = {'uniquely_mapped': uniquely_mapped,
                                            'multi_mapped': multi_mapped, 'mapped_sum': mapped_sum}
                        else:
                            reads[sj_id]['uniquely_mapped'] += uniquely_mapped
                            reads[sj_id]['multi_mapped'] += multi_mapped
                            reads[sj_id]['mapped_sum'] += mapped_sum

                        if reads[sj_id]['mapped_sum'] >= 20:
                            uniquely_mapped = reads[sj_id]['uniquely_mapped']
                            multi_mapped = reads[sj_id]['multi_mapped']
                            mapped_sum = reads[sj_id]['mapped_sum']
                            good_reads[sj_id] = {'uniquely_mapped': uniquely_mapped,
                                                 'multi_mapped': multi_mapped, 'mapped_sum': mapped_sum}
                    else:
                        good_reads[sj_id]['uniquely_mapped'] += uniquely_mapped
                        good_reads[sj_id]['multi_mapped'] += multi_mapped
                        good_reads[sj_id]['mapped_sum'] += mapped_sum

                else:  # mapped_sum >= 20
                    if sj_id not in good_reads.keys():
                        good_reads[sj_id] = {'uniquely_mapped': uniquely_mapped,
                                             'multi_mapped': multi_mapped, 'mapped_sum': mapped_sum}

                    else:  # already in good_reads
                        good_reads[sj_id]['uniquely_mapped'] += uniquely_mapped
                        good_reads[sj_id]['multi_mapped'] += multi_mapped
                        good_reads[sj_id]['mapped_sum'] += mapped_sum

            n += 1
# print(reads)
# print(good_reads)
print(f'Well suported introns number: {len(good_reads)}')

# Saving good SJ in csv file
df = pd.DataFrame(good_reads).transpose()
# print(df)
df.to_csv(f'{args.out_prefixName}.csv')
