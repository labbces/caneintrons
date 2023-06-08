# imports
import csv
import pandas as pd
import re
import argparse
import glob
from Bio import SeqIO

# Adding variables from terminal 
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta_file", type=str, required=True, help="All introns imeter type fasta file. Ex: data/introns/Souza_introns_IMEter/Souza_All_introns_imeter_intron.fa")
parser.add_argument("-s", "--STAR_path", type=str, required=True, help="Path to STAR output files. Ex: data/RNAseq_output/IQ_USP/Souza.sc.mlc.cns.sgl.utg.scga7.importdb.reformat_*_SJ.out.tab")
parser.add_argument("-n", "--filename", type=str, required=True, help="(Path to) output filename. Ex: SP803280_V2 - That would be saved as SP803280_V2.STAR_analysis.csv")
args = parser.parse_args()

fasta_file = args.fasta_file
star_files = args.STAR_path
filename = args.filename

# Dict to save fasta info
fasta_dict = {}

# Saving fasta file start and end coordinates per chromosome into a dict 
for record in SeqIO.parse(fasta_file, "fasta"):
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

# Final dict
STAR_dict = {}

# Going through STAR output TSV file
for starSRRfile in glob.glob(star_files):
    print(f'Current analyzing {starSRRfile}')
    with open(starSRRfile, 'r') as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter='\t')
        for row in tsv_reader:
            # from tsv getting chromosome, start, end uniquely and multimapped info
            chromosome, start, end, strand, _, annotated, uniquely_mapped, multi_mapped, _ = row
            uniquely_mapped = int(uniquely_mapped)
            if chromosome in fasta_dict.keys():
                #print(row)
                expression = int(uniquely_mapped) + int(multi_mapped)
                # Comparing STAR output start and end +/- 5 to position declared in the gtf from fasta file
                for s in range(int(start) - 5, int(start) + 6):
                    if s in fasta_dict[chromosome].keys():
                        for e in range(int(end) - 5, int(end) + 6):
                            if e in fasta_dict[chromosome][s].keys():
                                header = fasta_dict[chromosome][s][e]
                                # Filling final dict up
                                if header not in STAR_dict.keys():
                                    STAR_dict[header] = [expression, uniquely_mapped]
                                else:
                                    STAR_dict[header][0] += expression
                                    STAR_dict[header][1] += uniquely_mapped

print('STAR output analysis finished')
df = pd.DataFrame.from_dict(STAR_dict, orient='index', columns=['expression', 'uniquely_mapped'])
df = df.reset_index().rename(columns={'index': 'ID'})
# df = pd.DataFrame(list(STAR_dict.items()), columns=['ID', 'STAR_expression'])
df.to_csv(f'{filename}.STAR_analysis.csv')
print(f'{filename}.STAR_analysis.csv generated')        
