
import csv
import numpy as np
import matplotlib.pyplot as plt
# reads = [1, 2, 3 , 4, 5, 6,  10000, 10002, 10003, 4, 5, 6, 7, 1, 2, 3, 4, 5, 5, 6]


with open('/home/bia/sugarcane_introns_local/RNAseq_analysis/2output_CTBE.txt', 'r') as f:
    csv_string = f.read().strip()  # Lê o conteúdo do arquivo como uma string
    reads1 = csv_string.split(',')

reads = []
for element in reads1:
    reads.append(int(element))
print(min(reads))


media = np.mean(reads)
std_dev = np.std(reads)
median = np.median(reads)
var = np.var(reads)
quartis = np.percentile(reads, 25)
quartis2 = np.percentile(reads, 75)


print('---> SP80-3280-CTBE')
print(
    f'Média: {media}\nDesvio padrão: {std_dev}\nMediana: {median}\nVariância: {var}\nQuarts: (25) {quartis}\t(75)\t{quartis2}')


# print(reads)
reads = np.array(reads)
#reads = np.log10(reads)

with open('/home/bia/sugarcane_introns_local/RNAseq_analysis/2output_IQusp.txt', 'r') as f:
    csv_string = f.read().strip()  # Lê o conteúdo do arquivo como uma string
    reads1 = csv_string.split(',')

reads_IQ = []

for element in reads1:
    reads_IQ.append(int(element))
print(min(reads_IQ))


media = np.mean(reads_IQ)
std_dev = np.std(reads_IQ)
median = np.median(reads_IQ)
var = np.var(reads_IQ)
quartis = np.percentile(reads_IQ, 25)
quartis2 = np.percentile(reads_IQ, 75)

print('\n---> SP80-3280-IQ/USP')
print(
    f'Média: {media}\nDesvio padrão: {std_dev}\nMediana: {median}\nVariância: {var}\nQuarts: (25) {quartis}\t(75)\t{quartis2}')

# print(reads)
reads_IQ = np.array(reads_IQ)
#reads_IQ = np.log10(reads_IQ)


'''
plt.figure(figsize=(8, 6))
plt.boxplot(reads_IQ, labels=[''], positions=[3])
plt.yscale('log')
plt.ylabel('Quantidade de Splice Junctions')
plt.show()
'''

fig, ax = plt.subplots()
bp = ax.boxplot([reads_IQ, reads])
ax.set_xticklabels(['SP80-3280-IQ/USP', 'SP80-3280-CTBE'])
ax.set_yscale('log')
ax.set_ylabel('Quantidade de Leituras Mapeadas')
plt.show()
