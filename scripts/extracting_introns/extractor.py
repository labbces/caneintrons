# add argparse para os caminhos dos arquivos.

from Bio import SeqIO
import pyranges as pr
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--output", type=str, required=True)
args = parser.parse_args()


# Getting gff file
# gff = pr.read_gtf('/home/bia/sugarcane_introns/gtf_teste_cana.gtf')
gff = pr.read_gtf(
    '/home/bia/sugarcane_introns/3.Converting_to_gtf/from_trimmed_gff3\
/SP803280_171-181.gff3.gtf')

# GFF basic operations
# print(gff)
# print(gff.describe("transcript_id"))
# print(gff.features.introns(by="transcript"))

# Introns features pyranges and its basic operations
introns = gff.features.introns(by="transcript")
# print(introns)
# print(introns.Start)
# print(introns.End)

# Getting only introns that follows: 150 > intron size > 100
introns = (introns[(introns.End - introns.Start) > 99])
introns = (introns[(introns.End - introns.Start) < 151])
print(introns)

resume: dict = {}  # confirmar se isso está certo
for line_general in introns:
    for line in (line_general[1]['Chromosome']):
        id = line
        if id not in resume.keys():
            resume[id] = {'starts': [], 'ends': []}
    for line in (line_general[1]['Start']):
        resume[id]['starts'].append(line)
    for line in (line_general[1]['End']):
        resume[id]['ends'].append(line)
print(resume)
# print(resume.keys())

filename = args.output
with open('/home/bia/sugarcane_introns/GCA_002018215.1/ncbi_dataset/data/\
GCA_002018215.1/GCA_002018215.1_CTBE_SP803280_v1.0_genomic.fna') as handle:
    input_contigs = SeqIO.parse(handle, 'fasta')
    interest_contigs = (record for record in input_contigs if
                        record.id.split('.')[0] in resume.keys())
    print(input_contigs)
    print(interest_contigs)

# seqIO.index_db

    with open(f'/home/bia/sugarcane_introns/4_Extracting_introns/\
    {filename}.fa', 'a') as gene_file:
        SeqIO.write(interest_contigs, gene_file, "fasta")

    with open(f'/home/bia/sugarcane_introns/4_Extracting_introns/\
{filename}.fa', 'r') as gene_file:
        for record in SeqIO.parse(gene_file, "fasta"):
            seq_id = record.id.split('.')[0]
            # print(seq_id)
            seq_full = record.seq
            # colocar a posiçã original (sem a add para pegar as bordas dos introns)
            k = 0
            for start in resume[seq_id]['starts']:
                k += 1
                index = resume[seq_id]['starts'].index(start)
                end = (resume[seq_id]["ends"][index]) + 50
                start = start - 50
                intron_seq = seq_full[start:end]
                id_full = f'>{seq_id}.{k}'
                '''print(f'{id_full} ##### {intron_seq}')
                print(len(intron_seq) - 30)'''
                with open(f'/home/bia/sugarcane_introns/4_Extracting_introns/\
{filename}_intron.fa', 'a') as intron_file:
                    intron_file.write(f'{id_full}\n')
                    intron_file.write(f'{intron_seq}\n')
