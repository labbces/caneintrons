import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--chromosome", type=str, required=True,
                    help='')
args = parser.parse_args()


file = "data/Human_genome/knownAlt_hg18.txt"

chrs = []
with open(file, 'r') as infile:
    for line in infile:
        chr = line.split('\t')[1]
        if chr == args.chromosome:
            print(line.replace('\n', ''))
