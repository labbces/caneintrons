import logomaker
from Bio import SeqIO
import matplotlib.pyplot as plt

# Ler as sequências do arquivo FASTA
sequences = [str(record.seq) for record in SeqIO.parse("BothGenomes_spliceator_acceptor2.nr100.fa", "fasta")]

# Criar um objeto Logo
logo_data = logomaker.alignment_to_matrix(sequences, to_type="information")
logo = logomaker.Logo(logo_data, )
logo.ax.set_ylabel('Bits')
logo.ax.set_xlabel('Postions')
logo.ax.set_xticks(range(len(logo_data)))
#logo.ax.set_xticklabels('%+d'%x for x in [-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
logo.ax.set_xticklabels('%+d'%x for x in [-12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8])
logo.ax.axvline(11.5, color='k', linewidth=1, linestyle=':')
#logo.ax.axvline(9.5, color='k', linewidth=1, linestyle=':')
logo.ax.set_title('Acceptor splice site profile')
fig = plt.show()


