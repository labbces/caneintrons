from Bio import SeqIO
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--", type=str, required=True)
args = parser.parse_args()

fasta=args.fasta

contig_data = {}
for record in SeqIO.parse(fasta, "fasta"):
    contig_size = len(record.seq)
    contig_data[record.id] = contig_size

df = pd.DataFrame.from_dict(contig_data, orient="index")
df.to_csv("intron_size_data.csv")


