# imports
import pandas as pd
import numpy as np 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--spliceator_analysis", type=str, required=True, help="Path to spliceator analysis output file")
parser.add_argument("-d", "--splice2deep_analysis", type=str, required=True, help="Path to splice2deep analysis output file")
parser.add_argument("-s", "--star_analysis", type=str, required=True, help="Path to STAR analysis output file")
parser.add_argument("-a", "--IMEter_arabidopsis", type=str, required=True, help="Path to IMEter analysis output file - Arabidopsis thaliana model")
parser.add_argument("-o", "--IMEter_oriza", type=str, required=True, help="Path to IMEter analysis output file - Oriza model")
parser.add_argument("-b", "--IMEter_sorghum", type=str, required=True, help="Path to IMEter analysis output file - Sorghum model")
parser.add_argument("-f", "--filename", type=str, required=True, help="(Path to) output filename. Ex: SP803280 - that would be saved as SP803280.final_table.csv")
args = parser.parse_args()

spliceator_analysis = args.spliceator_analysis
splice2deep_analysis = args.splice2deep_analysis
star_analysis = args.star_analysis
IMEter_arabidopsis = args.IMEter_arabidopsis
IMEter_oriza = args.IMEter_oriza
IMEter_sorghum = args.IMEter_sorghum
filename = args.filename

df1 = pd.read_csv(spliceator_analysis) # Spliceator
df2 = pd.read_csv(splice2deep_analysis) # Splice2deep
df3 = pd.read_csv(IMEter_arabidopsis,sep='\t') # IMEter A thaliana
df4 = pd.read_csv(IMEter_oriza,sep='\t') # IMEter O sativa
df5 = pd.read_csv(IMEter_sorghum,sep='\t') # IMEter S bicolor
df6 = pd.read_csv(star_analysis, index_col=0) # STAR


'''
df1 = pd.read_csv('SP803280_D_.spliceator_analysis.csv') # Spliceator
df2 = pd.read_csv('Splice2Deep_analysis/analysis/CTBE_splice2deep_analysis.csv') # Splice2deep
df3 = pd.read_csv('/home/bia/sugarcane_introns_local/IMEter/IME_mycopy/results/CTBE_Athaliana_400_400.imeter',sep='\t') # IMEter A thaliana
df4 = pd.read_csv('/home/bia/sugarcane_introns_local/IMEter/IME_mycopy/results/CTBE_Osativa_400_400.imeter',sep='\t') # IMEter O sativa
df5 = pd.read_csv('/home/bia/sugarcane_introns_local/IMEter/IME_mycopy/results/CTBE_Sbicolor_400_400.imeter',sep='\t') # IMEter S bicolor
df6 = pd.read_csv('/home/bia/sugarcane_introns_local/RNAseq_analysis/analysis/CTBE.STARanalysis.csv') # STAR
'''

df1 = df1.drop("Unnamed: 0", axis=1)
columns = ['ID', 'Splice2Deep_Athaliana_Donor', 'Splice2Deep_Athaliana_Acceptor', 'Splice2Deep_Oriza_Donor', 'Splice2Deep_Oriza_Acceptor']
df2.columns = columns
columns = ['ID', 'IMEter_Athaliana_V1', 'IMEter_Athaliana_V2']
df3.columns = columns
columns = ['ID', 'IMEter_Osativa_V1', 'IMEter_Osativa_V2']
df4.columns = columns
columns = ['ID', 'IMEter_Sbicolor_V1', 'IMEter_Sbicolor_V2']
df5.columns = columns
df6['ID'] = df6['ID'].str.rstrip('\n')

'''
print(df1)
print(df2)
print(df3)
print(df4)
print(df5)
print(df6)
'''


df = df1.merge(df2, on='ID', how='inner')
df = df.merge(df3, on='ID', how='inner')
df = df.merge(df4, on='ID', how='inner')
df = df.merge(df5, on='ID', how='inner')
df = df.merge(df6, on='ID', how='inner')


df = df.dropna()
df = df.drop_duplicates()
df = df.reset_index(drop=True)

# df = pd.read_csv('/home/bia/sugarcane_introns_local/Final_Tables/CTBE_final_table_firstVersion.csv') # if final table already exist add here
# print(df)


df = df.drop(df[df['expression'] < 65].index)

df['Spliceator_average'] = (df['Spliceator_Donor_Score'] + df['Spliceator_Acceptor_Score']) / 2
min_value = df['Spliceator_average'].min()
max_value = df['Spliceator_average'].max()

df['Spliceator_normalize'] = (df['Spliceator_average'] - min_value) / (max_value - min_value)

min_value = df['IMEter_Athaliana_V2'].min()
max_value = df['IMEter_Athaliana_V2'].max()
df['IMEter_Athaliana_V2_normalize'] = (df['IMEter_Athaliana_V2'] - min_value) / (max_value - min_value)

min_value = df['IMEter_Sbicolor_V2'].min()
max_value = df['IMEter_Sbicolor_V2'].max()
df['IMEter_Sbicolor_V2_normalize'] = (df['IMEter_Sbicolor_V2'] - min_value) / (max_value - min_value)


sum = df['expression'].sum()
crm = (df['expression'] / sum) * (10 ** 6)
df['STAR_normalize'] = crm.apply(lambda x: np.log10(x))


df['s'] = (df['IMEter_Sbicolor_V2_normalize'] * 0.325) + (df['STAR_normalize'] * 0.2) + (df['Spliceator_normalize'] * 0.325) + (df['IMEter_Athaliana_V2_normalize'] * 0.12) + (1*0.03)

df = df.sort_values('s', ascending=False)
df = df.reset_index(drop=True)

# print(df)

df.to_csv(f'{filename}.final_table.csv')
# print(df)
