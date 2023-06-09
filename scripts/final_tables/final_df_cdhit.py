from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--cdhit_repseq", type=str, required=True, help="Path to cdhit representative sequences file")
parser.add_argument("-a", "--final_table1", type=str, required=True, help="Path to final_table from one genome")
parser.add_argument("-b", "--final_table2", type=str, required=True, help="Path to final_table from another genome")
parser.add_argument("-f", "--output_filename", type=str, required=True, help="(Path to) output filename. Ex: /home/usr/SP803280 - that will be saved as SP803280.cdhit_finalTable.csv")
args = parser.parse_args()

cdhit_representaive_output = args.cdhit_repseq
final_tableA = args.final_table1
final_tableB = args.final_table2
filename = args.output_filename

# cdhit_representaive_output = '/home/bia/sugarcane_introns_local/data/introns/cd-hit_output/BothGenomes_All_introns__300nr100'

representative_headers = []
for record in SeqIO.parse(cdhit_representaive_output, "fasta"):
    header = record.id
    print(header)
    representative_headers.append(header)

#df1 = pd.read_csv('/home/bia/sugarcane_introns_local/CTBE.final_table.csv', index_col=0)
#df2 = pd.read_csv('/home/bia/sugarcane_introns_local/Souza.final_table.csv', index_col=0)
df1 = pd.read_csv(final_tableA, index_col=0)
df2 = pd.read_csv(final_tableB, index_col=0)

df = pd.concat([df1, df2], ignore_index=True)
df = df[df['ID'].isin(representative_headers)]

#df1 = df1[df1['ID'].isin(representative_headers)]
#df2 = df2[df2['ID'].isin(representative_headers)]
#df = pd.concat([df1, df2], ignore_index=True)

df = df.sort_values('s', ascending=False)
df = df.reset_index(drop=True)

# print(df)
df.to_csv(f'{filename}.cdhit_finalTable.csv')
