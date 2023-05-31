
# SP = spliceator
# S2D = Splice2Deep
# O...= open ...

# imports
import argparse
import re
from Bio import SeqIO 
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--donor_ss_file", type=str, required=True, help="Path to file with donor splice sites")
parser.add_argument("-f", "--acceptor_ss_file", type=str, required=True, help="Path to file with acceptor splice sites")
parser.add_argument("-d", "--donor_s2d_file", type=str, required=True, help="Path to donor Splice2Deep output files")
parser.add_argument("-a", "--acceptor_s2d_file", type=str, required=True, help="Path to acceptor Splice2Deep output files")
parser.add_argument("-p", "--spliceator_introns_donor", type=str, required=True, help="Path to file containing all good donors predicted by spliceator")
parser.add_argument("-t", "--spliceator_introns_acceptor", type=str, required=True, help="Path to file containing all good donors predicted by spliceator")
parser.add_argument("-i", "--all_introns", type=str, required=True, help="Path to full length sugarcane introns + borders")

args = parser.parse_args()

donor_ss_file = args.donor_ss_file
acceptor_ss_file = args.acceptor_ss_file
donor_s2d_file = args.donor_s2d_file
acceptor_s2d_file = args.acceptor_s2d_file
Dspliceator_introns = args.spliceator_introns_donor
Aspliceator_introns = args.spliceator_introns_acceptor
all_introns = args.all_introns


with open(donor_ss_file, 'r') as Odonor_ss_file, open(donor_s2d_file, 'r') as Odonor_s2d_file, open(acceptor_s2d_file, 'r') as Oacceptor_s2d_file:

    data = {}

    # Usar o SeqIO para ler o arquivo FASTA e iterar pelas linhas dos três arquivos simultaneamente
    for seq_record, linha2, linha3 in zip(SeqIO.parse(Odonor_ss_file, 'fasta'), Odonor_s2d_file, Oacceptor_s2d_file):

        # Extrair o ID da sequência
        seq_id = seq_record.id
        data[seq_id] = [linha2.strip(), linha3.strip()]

#print(data)

df = pd.DataFrame.from_dict(data, orient='index', columns=['Donor', 'Acceptor'])
#print(df)
df = df.loc[(df['Donor'] == '1') & (df['Acceptor'] == '1')]
#print(df)

# donor_ss_file
#fasta = all_introns
#inx = fasta+'.genome_inx'
#fasta_index = SeqIO.index_db(inx, fasta, 'fasta')

# Donor spliceator
fasta = Dspliceator_introns
#fasta = '/home/bia/sugarcane_introns_local/spliceator_anaysis/BothGenomes_spliceator_donor2.fa'
inx = fasta+'.genome_inx'
fasta_index_spD = SeqIO.index_db(inx, fasta, 'fasta')

# Acceptor splcieator
fasta = Aspliceator_introns
#fasta = '/home/bia/sugarcane_introns_local/spliceator_anaysis/BothGenomes_spliceator_acceptor2.fa'
inx = fasta+'.genome_inx'
fasta_index_spA = SeqIO.index_db(inx, fasta, 'fasta')

# Splice2Deep donor
#fasta = donor_ss_file
#inx = fasta+'.genome_inx'
#fasta_index_s2dD = SeqIO.index_db(inx, fasta, 'fasta')

# Splice2Deep acceptor
#fasta = acceptor_ss_file
#inx = fasta+'.genome_inx'
#fasta_index_s2dA = SeqIO.index_db(inx, fasta, 'fasta')



common_ids = 0
good_introns = 0
m = re.match(r'(.*_[0-9]*-[0-9]*).*[acceptor/donor]\_intron.fa', donor_ss_file)
fasta_output = m.group(1)+'_splice2deep_goodIntrons.fasta'
fasta_outputA = m.group(1)+'_splice2deep_goodAcceptors.fasta'
fasta_outputD = m.group(1)+'_splice2deep_goodDonors.fasta'

with open(fasta_output, 'w') as output, open(fasta_outputA,'w') as fasta_outputA, open(fasta_outputD,'w') as fasta_outputD:
    for index, row in df.iterrows():
        good_introns += 1

        # seq_full = fasta_index[index].seq


        #print(index)
        #print(f'>{index}')
        #print(seq_full)

        # donor_seq = fasta_index_s2dD[index].seq
        # acceptor_seq = fasta_index_s2dA[index].seq


        # output.write(f'>{index}\n')
        # output.write(f'{str(seq_full)}\n')

        #fasta_outputA.write(f'>{index}\n')
        #fasta_outputA.write(f'{str(acceptor_seq)}\n')

        # fasta_outputD.write(f'>{index}\n')
        # fasta_outputD.write(f'{str(donor_seq)}\n')

        if index in fasta_index_spD and index in fasta_index_spA:
        # if index in fasta_index_spD:
            common_ids += 1

print(f'Number of shared ids - spliceator and splice2deep: {common_ids}')
print(f'Number of good introns - splice2deep: {good_introns}')

#  python3 -u "/home/bia/sugarcane_introns_local/merge_splice2deep.py" -f /home/bia/sugarcane_introns_local/SP803280_D_151-161.gtfacceptor_intron.fa -d /home/bia/sugarcane_introns_local/data/introns/CTBE_introns_s2d/splicedeep_DoSS_output_SP803280_D_151-161.gtfdonor_intron.fa.splice2deep -a /home/bia/sugarcane_introns_local/data/introns/CTBE_introns_s2d/splicedeep_AcSS_output_SP803280_D_151-161.gtfacceptor_intron.fa.splice2deep -s /home/bia/sugarcane_introns_local/SP803280_D_151-161.gtfdonor_intron.fa -p /home/bia/sugarcane_introns_local/spliceator_anaysis/BothGenomes_spliceator_donor2.fa -t /home/bia/sugarcane_introns_local/spliceator_anaysis/BothGenomes_spliceator_acceptor2.fa -i /home/bia/sugarcane_introns_local/data/introns/CTBE_introns/SP803280_D_151-161.gtf_intron.fa
