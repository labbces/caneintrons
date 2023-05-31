'''
Script to analyze spliceator output files.

Now the data is splited into 2 types of files:
- One has plus and bad generated neg strands
- Other only with good generated plus strand

Objetcive: Generate a csv file with good introns ID, Donor Score and Acceptor score
predicted by spliceator tool 
'''

# Imports
import pandas as pd
import re
import argparse
import glob
import os

# Unisng argparse to simplify files loading to script
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--spliceator_plus_output_files", type=str, required=True, help="Path too directory with all SRR star output files")
parser.add_argument("-n", "--spliceator_neg_output_files", type=str, required=True, help="Path too directory with all SRR star output files")
parser.add_argument("-i", "--ids", type=str, required=True, help="file with good introns IDs")
parser.add_argument("-f", "--output_filename", type=str, required=True, help="XXXXXXX.spliceatorAnalysis.csv where XXXXX id the -f string")
args = parser.parse_args()

# Creating variables from argparse
plus_outfile = args.spliceator_plus_output_files
neg_outfile =args.spliceator_neg_output_files
filename = args.output_filename
ids = args.ids

verify = []
with open(ids, 'r') as ids:
    for line in ids:
        if line.startswith('>'):
            verify.append(line.lstrip(">").rstrip("\n"))
# print(verify)

 
#for splcieator_file in glob.glob(spliceator_plus_output_files):
with open(plus_outfile, 'r') as plus, open(neg_outfile, 'r') as neg:    
    df_plus = pd.read_csv(plus, sep=';')
    # print(df_plus)
    df_plus_only = df_plus[~df_plus['ID'].str.contains('---')]
    #  print(df_plus_only)

    mask = df_plus_only['ID'].isin(verify)
    df_plus = df_plus_only[mask]
    #print(df_plus) 

    df_plus_donor = df_plus.loc[df_plus['SS_type'] == "Donor"]
    df_plus_acceptor = df_plus.loc[df_plus['SS_type'] == "Acceptor"]
    #print(df_plus_donor)
    #print(df_plus_acceptor)


    df_neg = pd.read_csv(neg, sep=';')
    #print(df_neg)
    mask = df_neg['ID'].isin(verify)
    df_neg = df_neg[mask]
    print(df_neg)

    df_neg_donor = df_neg.loc[df_neg['SS_type'] == "Donor"]
    df_neg_acceptor = df_neg.loc[df_neg['SS_type'] == "Acceptor"]

    df_joined_plus = pd.merge(df_plus_donor, df_plus_acceptor, on="ID", how="inner")
    df_joined_neg = pd.merge(df_neg_donor, df_neg_acceptor, on="ID", how="inner")

    colunas_desejadas = ['ID', 'Score_x', 'Score_y']
    merged_df_plus = df_joined_plus[colunas_desejadas]
    merged_df_neg = df_joined_neg[colunas_desejadas]

    #print(merged_df)
    merged_df_plus.columns = ['ID', 'Spliceator_Donor_Score', 'Spliceator_Acceptor_Score']
    merged_df_neg.columns = ['ID', 'Spliceator_Donor_Score', 'Spliceator_Acceptor_Score']

    df_concatenado = pd.concat([merged_df_plus, merged_df_neg], axis=0)
    print(df_concatenado)


    df_concatenado['Media'] = df_concatenado[['Spliceator_Donor_Score', 'Spliceator_Acceptor_Score']].mean(axis=1)

    # Ordenar o DataFrame pelo valor da coluna 'Media' em ordem decrescente
    df_concatenado = df_concatenado.sort_values('Media', ascending=False)

    # Criar uma visualização correta do DataFrame
    # Manter apenas as linhas únicas com maior valor na coluna 'Media' para cada ID
    df_unique = df_concatenado.drop_duplicates('ID', keep='first')

    # Excluir a coluna 'Media'
    df_unique = df_unique.drop('Media', axis=1)
    #print(df_unique)

    df_unique = df_unique.reset_index(drop=True)


df_unique.to_csv('SP803280_D_.spliceator_analysis.csv')
# merged_df.to_csv('TESTE2')

# /home/bia/sugarcane_introns_local/spliceator_anaysis/spliceatorOutput/CTBE/
