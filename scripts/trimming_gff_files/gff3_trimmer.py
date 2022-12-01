import statistics as st
import math
import argparse

# /home/bia/sugarcane_introns/old_nov30/CTBE_SP803280_genome_annotation.EVM19022016.fixIDS.gff3


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--gtf_file", type=str, required=True)
parser.add_argument("-o", "--output_base_name",
                    type=str, required=True)  # variedade
parser.add_argument("-p", "--path_to_save", type=str, default='./')
args = parser.parse_args()

gtf = args.gtf_file
base_name = args.output_base_name
path = args.path_to_save

if path[-1] != "/":
    path = path+'/'

# Getting contigs frequencies
contigs_freq = {}
contigs_line = {}

with open(gtf, 'r') as gtf:
    for line in gtf:
        contig = line.split('\t')[0]
        try:
            contigs_freq[contig] += 1
        except:
            contigs_freq[contig] = 1

        try:
            contigs_line[contig].append(line)
        except:
            contigs_line[contig] = [line]
# print(contigs_freq)
# print(contigs_line)


freq = list(contigs_freq.values())
# print(freq)


# Getting statistics
mean = st.mean(freq)
print(f'Média: {mean}')

mediana = st.median(freq)
print(f'Mediana: {mediana}')

max = max(freq)
print(f'Máximo {max}')

min = min(freq)
print(f'Mínimo: {min}')

amplitude = max - min
# print(amplitude)

n = len(freq)
print(f'Quantidade de genes: {n}')

# Getting intervals
k = 1 + 3.332 * math.log(n, 10)
print(f'Quantidade de classes: {k}')

amplitude_do_intervalo = amplitude/k
print(f'Amplitude do intervalo: {amplitude_do_intervalo}')

# Separating classes
contador = 0
todos_os_contigs = []
while contador < int(k) + 1:
    inferior_limit = min + contador*amplitude_do_intervalo
    upper_limit = inferior_limit + amplitude_do_intervalo

    classe = [contig for contig in contigs_freq.keys() if contigs_freq[contig]
              >= inferior_limit and contigs_freq[contig] < upper_limit]
    # print(len(classe))

    filename = f'{base_name}_{int(inferior_limit)}-{int(upper_limit)}.gtf'

    if len(classe) > 0:
        for contig in classe:
            if contig in todos_os_contigs:
                print(
                    f'ERRO: código do contig aparece em\
 mais de uma classe {contig}')
                quit()
            else:
                todos_os_contigs.append(contig)
                for line in contigs_line[contig]:
                    with open(f'{path}{filename}', 'a') as newfile:
                        newfile.write(line)
    contador += 1
