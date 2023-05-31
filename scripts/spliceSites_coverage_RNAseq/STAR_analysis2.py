'''import pandas as pd
import re
import argparse
import glob
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--good_introns_file", type=str, required=True, help="Path to good introns file")
parser.add_argument("-s", "--star_output_files", type=str, required=True, help="Path too directory with all SRR star output files")
parser.add_argument("-f", "--output_filename", type=str, required=True, help="XXXXXXX.STARanalysis.csv where XXXXX id the -f string")
args = parser.parse_args()

good_introns_file_ = args.good_introns_file
star_output_path = args.star_output_files
filename = args.output_filename



starDF = {'ID': [], 'star_expression': []}
info_dict = {}

for record in SeqIO.parse(good_introns_file_, "fasta"):
    header = record.id
    header = header.lstrip(">")
    m = re.match('([A-Za-z_.a0-9._]*)\.([0-9]*)-([0-9]*)-([\+\-])-(.*)', header)
    chr = m.group(1)
    start = int(m.group(2)) + 1
    end = int(m.group(3))
    # print(f"Header: {header}\nchr: {chr}\nstart: {start}\n end: {end}")
    if chr not in info_dict.keys():
        info_dict[chr] = {header: {"start": start, 'end': end}}
    else:
        info_dict[chr][header] = {"start": start, 'end': end}
print(f'Genome processed...\n')

# print(info_dict)

for starSRRfile in glob.glob(star_output_path):
    print(f'Processing {starSRRfile}...')
    df = pd.read_csv(starSRRfile, delimiter='\t', header=None)
    columns = ['chromosome', 'start', 'end', 'strand', 'splice_site', 'annotated', 'uniquely_mapped', 'multi_mapped', 'overhang']
    df.columns = columns

    b = 0
    for chr in info_dict.keys():
        b += 1
        #print(f'\tProcessing chromosome {chr}...')
        percentage = b/len(info_dict)*100
        print(f'\t{b}/{len(info_dict)} - {percentage:.2f}%')
        df_chr = df.loc[df['chromosome'] == chr]
        l = 0
        for header in info_dict[chr].keys():
            l += 1
            # print(f'\t\t{l}/{len(info_dict[chr])} - {(l/len(info_dict[chr]))*100:.2f}%')
            start = info_dict[chr][header]['start']
            end = info_dict[chr][header]['end']
            df_good_introns = df_chr.loc[(df_chr['start'] >= start - 5) & (df_chr['start'] <= start + 5) &
                                         (df_chr['end'] >= end - 5) & (df_chr['end'] <= end + 5)]
            num_lines = df_good_introns.shape[0]
            if num_lines > 0:
                x = 0
                while x < num_lines:
                    expression = df_good_introns.loc[df_good_introns.index[x], 'uniquely_mapped'] + df_good_introns.loc[df_good_introns.index[x], 'multi_mapped']
                    if header not in starDF['ID']:
                        starDF['ID'].append(header)
                        starDF['star_expression'].append(expression)
                    else:
                        pos = starDF['ID'].index(header)
                        starDF['star_expression'][pos] += expression
                    x += 1
        df = df.drop(df_chr.index)

starDF = pd.DataFrame(starDF)
filename = f'{filename}.STARanalysis.csv'
starDF.to_csv(filename, index=False)





import re
# Carregar o arquivo fasta em um dicionário
fasta_dict = {}
with open('data/introns/Souza_introns_IMEter/Souza_All_introns_imeter_intron.fa', 'r') as fasta_file:
    header = ''
    sequence = ''
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if header and sequence:
                fasta_dict[header] = sequence
            header = line[1:]
            sequence = ''
        else:
            sequence += line
    # Adicionar a última sequência do arquivo
    if header and sequence:
        fasta_dict[header] = sequence

# Criar dicionário para armazenar a soma dos valores
expression_dict = {}

'''


import csv
import pandas as pd
import re
import argparse
import glob
import os
from Bio import SeqIO

fasta_dict = {}

for record in SeqIO.parse('data/introns/Souza_introns_IMEter/Souza_All_introns_imeter_intron.fa', "fasta"):
    header = record.id
    header = header.lstrip(">")
    m = re.match('([A-Za-z_.a0-9._]*)\.([0-9]*)-([0-9]*)-([\+\-])-(.*)', header)
    chr = m.group(1)
    start = int(m.group(2)) 
    end = int(m.group(3))
    if chr not in fasta_dict:
        fasta_dict[chr] = {start: {end: header}}
    else:
        fasta_dict[chr][start] = {end: header}

print(f'Fasta Processed...')
#print(fasta_dict)

STAR_dict = {}
# Percorrer o arquivo tsv

for starSRRfile in glob.glob('data/RNAseq_output/IQ_USP/Souza.sc.mlc.cns.sgl.utg.scga7.importdb.reformat_*_SJ.out.tab'):
    print(f'Current analyzing {starSRRfile}')
    with open(starSRRfile, 'r') as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter='\t')
        for row in tsv_reader:
            chromosome, start, end, strand, _, annotated, uniquely_mapped, multi_mapped, _ = row
            if chromosome in fasta_dict.keys():
                #print(row)
                expression = int(uniquely_mapped) + int(multi_mapped)
                for s in range(int(start) - 5, int(start) + 6):
                    #print(f's {s}')
                    if s in fasta_dict[chromosome].keys():
                        for e in range(int(end) - 5, int(end) + 6):
                            if e in fasta_dict[chromosome][s].keys():
                                header = fasta_dict[chromosome][s][e]
                                #print(f'HEADER: {header}')
                                #print(f'{s} {e}')

                                if header not in STAR_dict.keys():
                                    STAR_dict[header] = expression
                                else:
                                    STAR_dict[header] += expression

    # print(STAR_dict)

df = pd.DataFrame(list(STAR_dict.items()), columns=['ID', 'STAR_expression'])
df.to_csv('Souza_STAR_analysis.csv')
        

