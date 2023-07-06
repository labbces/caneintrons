from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--chromosome", type=str, required=True,
                    help='Chromosome name - to make the analysis faster')
parser.add_argument("-g", "--gtf", type=str,
                    required=True, help='genome gtf file')
args = parser.parse_args()

# create a dictionary of dsictionaries using defaultdict
# outer dictionary key is gene_id
# inner dictionary key is transcript_id

data = defaultdict(dict)

# create a default dic with two inner dictionaries
# outer key is gene_id
# inner key is exon position
# inner value is another dictionary with transcript_id as key
exonPos = defaultdict(lambda: defaultdict(dict))
exonPos2 = defaultdict(lambda: defaultdict(dict))

file = args.gtf
# file = "/home/bia/sugarcane_introns_local/data/GTF/Homo_sapiens.NCBI36.54.gtf"
# file = "/home/bia/sugarcane_introns_local/data/GTF/hg38.ensGene.gtf"
# open a gtf file
with open(file, "r") as infile:
    # read the file line by line
    for line in infile:
        # remove the newline character
        line = line.strip()
        #skip lines starting with #
        if line.startswith("#"):
            continue
        # split the line into a list
        fields = line.split("\t")
        # check if the genome version is correct
        # if fields[1].upper().strip() == "TAIR10":
        # check if field[2] is a trranscript
        if fields[2] == "exon":
            strand = fields[6]
            ninethfiled = fields[8].split(";")
            for tag in ninethfiled:
                # rermove leading spaces
                tag = tag.strip()
                if tag.startswith("transcript_id"):
                    tid = tag.split()[1]
                elif tag.startswith("gene_id"):
                    gid = tag.split()[1]
            if gid not in data:
                data[gid] = {}
            if tid not in data[gid]:
                data[gid][tid] = 0
            # count the exons
            exonPositions = str(
                fields[0]+':'+fields[3]+"-"+fields[4]+":"+strand)
            exonPositions2 = [fields[0], int(
                fields[3]), int(fields[4]), strand]
            exonPos[gid][exonPositions][tid] = 0
            # check exon overlapping
            for id, exon in exonPos2[gid][tid].items():
                if int(fields[3]) <= exon[2]:  # overlap
                    if len(exonPositions2) == 4:
                        exonPositions2.append('overlap')
                    if len(exonPos2[gid][tid][id]) == 4:
                        exonPos2[gid][tid][id].append('overlap')
            exonPos2[gid][tid][exonPositions] = exonPositions2
            data[gid][tid] += 1

a = []
# print('geneID', 'transcriptID', 'coordinate', 'transcript_Amount', 'exons_Amount', sep="\t")
for gid in data.keys():
    if len(data[gid]) >= 4:
        for tid in data[gid].keys():
            if data[gid][tid] >= 3:
                for exon, pos in exonPos2[gid][tid].items():
                    if len(pos) == 4:
                        if len(exonPos[gid][exon]) == len(data[gid]):
                            #chr = exon.split(':')[0]
                            # if chr.strip() == args.chromosome:
                            #print(gid, tid, exon, len(data[gid]), data[gid][tid], sep="\t")
                            type_id = str(tid).replace('"', '')
                            print(f'CS_{type_id}', f'chr{exon}', sep='\t')
