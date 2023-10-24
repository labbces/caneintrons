import pandas as pd
import numpy as np

# unique sequences
cdhitFASTA = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/newIntrons/fastas/Both_newIntrons_IMEter_intron_nr100.fa', header=None)
cdhitFASTA = cdhitFASTA[cdhitFASTA.iloc[:, 0].str.startswith('>')]
cdhitFASTA.columns = ['ID']
# preditores
df1 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/spliceator_anaysis/newIntrons/Both_newIntrons_spliceator_analysis.csv')  # Spliceator
df1.columns = ['ID', 'donor_sequence',
               'acceptor_sequence', 'strand', 'Spliceator_average']
df2 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/Splice2Deep/newIntrons/Both_newIntrons_s2d_analysis.csv')  # Splice2deep

# IMEter
# IMEter A thaliana


def adicionar_maior_que(text):
    return ">" + text


df3 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/IMEter/IME_mycopy/results/Both_newIntrons_Athaliana_400_400.imeter', sep='\t', header=None)
df3.columns = ['ID', 'IMEter_Athaliana_V1', 'IMEter_Athaliana_V2']
df3['ID'] = df3['ID'].apply(adicionar_maior_que)
# IMEter O sativa
df4 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/IMEter/IME_mycopy/results/Both_newIntrons_Osativa_400_400.imeter', sep='\t', header=None)
df4.columns = ['ID', 'IMEter_Osativa_V1', 'IMEter_Osativa_V2']
df4['ID'] = df4['ID'].apply(adicionar_maior_que)
# IMEter S bicolor
df5 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/IMEter/IME_mycopy/results/Both_newIntrons_Sbicolor_400_400.imeter', sep='\t', header=None)
df5.columns = ['ID', 'IMEter_Sbicolor_V1', 'IMEter_Sbicolor_V2']
df5['ID'] = df5['ID'].apply(adicionar_maior_que)
# STAR
df6 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/newIntrons/STAR/Both_NewIntrons_STAR_filtered.csv')  # STAR
df6['ID'] = '>' + df6['contig'] + '.' + \
    df6['start'].astype(str) + '-' + \
    (df6['end'] + 1).astype(str) + '-' + df6['strand']

#print(cdhitFASTA, df1, df2, df3, df4, df5, df6)
print(df6)

dataframes = [cdhitFASTA, df1, df2, df3, df4, df5, df6]
df = dataframes[0]
for dfs in dataframes[1:]:
    df = df.merge(dfs, on='ID', how='outer')
df = df.dropna()
print(df)
# Save the DataFrame to a CSV file without including the index

min_value = df['Spliceator_average'].min()
max_value = df['Spliceator_average'].max()
df['Spliceator_normalized'] = (
    df['Spliceator_average'] - min_value) / (max_value - min_value)

min_value = df['IMEter_Athaliana_V2'].min()
max_value = df['IMEter_Athaliana_V2'].max()
df['IMEter_Athaliana_V2_normalized'] = (
    df['IMEter_Athaliana_V2'] - min_value) / (max_value - min_value)

min_value = df['IMEter_Sbicolor_V2'].min()
max_value = df['IMEter_Sbicolor_V2'].max()
df['IMEter_Sbicolor_V2_normalized'] = (
    df['IMEter_Sbicolor_V2'] - min_value) / (max_value - min_value)


sum = df['mapped_sum'].sum()
print(sum)
cpm = (df['mapped_sum'] / sum) * (10 ** 6)
df['STAR_log10CPM'] = cpm.apply(lambda x: np.log10(x))
min_value = df['STAR_log10CPM'].min()
max_value = df['STAR_log10CPM'].max()
df['STAR_normalized'] = (df['STAR_log10CPM'] -
                         min_value) / (max_value - min_value)
print(df)
df['s'] = (df['IMEter_Sbicolor_V2_normalized'] * 0.325) + (df['STAR_normalized'] * 0.27) + \
    (df['Spliceator_normalized'] * 0.325) + \
    (df['IMEter_Athaliana_V2_normalized'] * 0.05) + (1*0.03)

df = df.sort_values('s', ascending=False)
df = df.reset_index(drop=True)
df.to_csv('output_file.csv', index=True)
