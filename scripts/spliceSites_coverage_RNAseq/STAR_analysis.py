

import pandas as pd
import re
import argparse
import glob
import os

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--good_introns_file", type=str, required=True, help="Path to good introns file")
parser.add_argument("-s", "--star_output_files", type=str, required=True, help="Path too directory with all SRR star output files")
parser.add_argument("-f", "--output_filename", type=str, required=True, help="XXXXXXX.STARanalysis.csv where XXXXX id the -f string")
args = parser.parse_args()

good_introns_file_ = args.good_introns_file
star_output_path = args.star_output_files
filename = args.output_filename



starDF = {'ID': [], 'star_expression': []}

# '/home/bia/sugarcane_introns_local/data/RNAseq_output/GCA_CTBE/GCA_002018215.1_CTBE_SP803280_v1.0_genomic_*_SJ.out.tab'
for starSRRfile in glob.glob(star_output_path):
    df = pd.read_csv(starSRRfile, delimiter='\t', header=None)
    columns = ['chromosome', 'start', 'end', 'strand', 'splice_site', 'annotated', 'uniquely_mapped', 'multi_mapped', 'overhang']
    df.columns = columns
    
    df_num = df.shape[0]
    print(f'Processing {starSRRfile}... {df_num} records\n')
    k = 0

    with open(good_introns_file_, 'r') as good_introns_file:
        for good_intron in good_introns_file:
            if good_intron.startswith('>'):
                k += 1
                print(k)
                print(good_intron)
                header = good_intron.lstrip(">")
                m = re.match('>([A-Za-z_.a0-9._]*)\.([0-9]*)-([0-9]*)-([\+\-])-(.*)', good_intron)
                chr = m.group(1)
                start = int(m.group(2)) + 1
                end = int(m.group(3))

            
                subdf = df.loc[(df['chromosome'] == chr) & 
                                (df['start'] >= start - 5) & (df['start'] <= start + 5) &
                                (df['end'] >= end - 5) & (df['end'] <= end + 5)]
                
                num_lines = subdf.shape[0]
                if num_lines > 0:
                    #print(subdf)
                    x = 0
                    while x < num_lines:
                        print(f'teste {x}')
                        expression = subdf.loc[subdf.index[x], 'uniquely_mapped'] + subdf.loc[subdf.index[x], 'multi_mapped']

                        if header not in starDF['ID']:
                            starDF['ID'].append(header)
                            starDF['star_expression'].append(expression)
                        else:
                            pos = starDF['ID'].index(header)
                            starDF['star_expression'][pos] += expression
                        x += 1
                df = df.drop(subdf.index)

starDF = pd.DataFrame(starDF)
filename = f'{filename}.STARanalysis.csv'
starDF.to_csv(filename, index=False)


#python3 -u "/home/bia/sugarcane_introns_local/STAR_analysis.py" -g /home/bia/sugarcane_introns_local/data/introns/CTBE_introns/CTBE_All_introns.fa -s '/home/bia/sugarcane_introns_local/data/RNAseq_output/GCA_CTBE/GCA_002018215.1_CTBE_SP803280_v1.0_genomic_*_SJ.out.tab' -f TESTE_CTBE
# is in