from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

cdhit_representaive_output = '/home/bia/sugarcane_introns_local/data/introns/cd-hit_output/BothGenomes_All_introns__300nr100'

representative_headers = []
for record in SeqIO.parse(cdhit_representaive_output, "fasta"):
    header = record.id
    representative_headers.append(header)

df1 = pd.read_csv('/home/bia/sugarcane_introns_local/CTBE.final_table.csv', index_col=0)
df2 = pd.read_csv('/home/bia/sugarcane_introns_local/Souza.final_table.csv', index_col=0)
df1 = df1[df1['ID'].isin(representative_headers)]
df2 = df2[df2['ID'].isin(representative_headers)]

df = pd.concat([df1, df2], ignore_index=True)
df = df.sort_values('s', ascending=False)
df = df.reset_index(drop=True)

print(df)
df.to_csv('filename.csv')
