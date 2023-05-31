
# SP = spliceator
# S2D = Splice2Deep
# O...= open ...

# imports
import argparse
import re
from Bio import SeqIO 
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--donor_ss_file", type=str, required=True, help="Path to file with donor splice sites")
parser.add_argument("-f", "--acceptor_ss_file", type=str, required=True, help="Path to file with acceptor splice sites")
parser.add_argument("-d", "--donor_s2d_file", type=str, required=True, help="Path to donor Splice2Deep output files")
parser.add_argument("-a", "--acceptor_s2d_file", type=str, required=True, help="Path to acceptor Splice2Deep output files")
parser.add_argument("-p", "--spliceator_introns_donor", type=str, required=True, help="Path to file containing all good donors predicted by spliceator")
parser.add_argument("-t", "--spliceator_introns_acceptor", type=str, required=True, help="Path to file containing all good donors predicted by spliceator")
parser.add_argument("-i", "--all_introns", type=str, required=True, help="Path to full length sugarcane introns + borders")
parser.add_argument("-b", "--orizaModel_Donor_s2d_Output", type=str, help="")
parser.add_argument("-c", "--orizaModel_Acceptor_s2d_Output", type=str, help="")
parser.add_argument("-n", "--filename", type=str, help="")


args = parser.parse_args()

# Defining variables
donor_ss_file = args.donor_ss_file
acceptor_ss_file = args.acceptor_ss_file
donor_s2d_file = args.donor_s2d_file
acceptor_s2d_file = args.acceptor_s2d_file
Dspliceator_introns = args.spliceator_introns_donor
Aspliceator_introns = args.spliceator_introns_acceptor
all_introns = args.all_introns
orizaDonor = args.orizaModel_Donor_s2d_Output
orizaAcceptor = args.orizaModel_Acceptor_s2d_Output

# Openning acceptor and donor splice sites (ss) and s2d output to correlate and saving it to correlate to spliceaator output 
with open(donor_ss_file, 'r') as Odonor_ss_file, open(donor_s2d_file, 'r') as Odonor_s2d_file, open(acceptor_s2d_file, 'r') as Oacceptor_s2d_file, open(orizaDonor, 'r') as OorizaDonor, open(orizaAcceptor, 'r') as OorizaAcceptor:

    data = {}

    # Usar o SeqIO para ler o arquivo FASTA e iterar pelas linhas dos três arquivos simultaneamente
    for seq_record, linha2, linha3, linha4, linha5 in zip(SeqIO.parse(Odonor_ss_file, 'fasta'), Odonor_s2d_file, Oacceptor_s2d_file, OorizaDonor, OorizaAcceptor):

        # Extrair o ID da sequência
        seq_id = seq_record.id
        data[seq_id] = [linha2.strip(), linha3.strip(), linha4.strip(), linha5.strip()]

# DataFrame with all s2d output data
dframe = pd.DataFrame.from_dict(data, orient='index', columns=['Donorat', 'Acceptorat', 'DonorOr', 'AcceptorOr'])



print(dframe.shape[0])
df = dframe.loc[((dframe['Donorat'].astype(int) == 1) & (dframe['Acceptorat'].astype(int) == 1)) | ((dframe['DonorOr'].astype(int) == 1) & (dframe['AcceptorOr'].astype(int) == 1))]

analysis_output = args.filename +'_splice2deep_analysis.csv'
df.to_csv(analysis_output)

good_introns_ar = df.loc[(df['Donorat'].astype(int) == 1) & (df['Acceptorat'].astype(int) == 1)]
good_introns_Oriza = df.loc[(df['DonorOr'].astype(int) == 1) & (df['AcceptorOr'].astype(int) == 1)]



good_introns_both = df.loc[(df['DonorOr'].astype(int) == 1) & (df['AcceptorOr'].astype(int) == 1) & (df['Donorat'].astype(int) == 1) & (df['Acceptorat'].astype(int) == 1)]

# print(df)


fastaAll = all_introns
inxAll = fastaAll +'.genome_inx'
fasta_index_All = SeqIO.index_db(inxAll, fastaAll, 'fasta')

# Donor spliceator
fasta = Dspliceator_introns
inx = fasta+'.genome_inx'
fasta_index_spD = SeqIO.index_db(inx, fasta, 'fasta')

# Acceptor splcieator
fasta = Aspliceator_introns
inx = fasta+'.genome_inx'
fasta_index_spA = SeqIO.index_db(inx, fasta, 'fasta')

common_ids = 0
common_ids_or = good_introns_Oriza.shape[0]
common_ids_At = good_introns_ar.shape[0]
common_ids_both = good_introns_both.shape[0]

good_introns_s2d = df.shape[0]


fasta_output = args.filename +'_splice2deep_goodIntrons.fasta'

#fasta_outputA = m.group(1)+'_splice2deep_goodAcceptors.fasta'
#fasta_outputD = m.group(1)+'_splice2deep_goodDonors.fasta'

good_headers = []

with open(fasta_output, 'w') as output:
    for index in df.index:
        if index in fasta_index_spD and index in fasta_index_spA:
            common_ids += 1
            seq_full = fasta_index_All[index].seq
            output.write(f'>{index}\n')
            output.write(f'{seq_full}\n')


print(f'Number of good introns according to Splice2deep (At or Oriza): {good_introns_s2d}')
print(f'Number of good introns according to Splice2deep (At and Oriza): {common_ids_both}')
print(f'Number of good introns according to Splice2deep (At):  {common_ids_At}')
print(f'Number of good introns according to Splice2deep (Oriza):  {common_ids_or}')
print(f'Number of shared ids - spliceator and splice2deep: {common_ids}')


# python3 -u "/home/bia/sugarcane_introns_local/merge_splice2deep_2.py" -s data/introns/CTBE_introns_s2d/SP803280_D_00-001.gtfdonor_intron.fa -f data/introns/CTBE_introns_s2d/SP803280_D_00-001.gtfacceptor_intron.fa -d Splice2Deep_analysis/output_s2d/CTBE/splicedeep_DoSS_output_SP803280_D_00-001.gtfdonor_at_intron.fa.splice2deep -a Splice2Deep_analysis/output_s2d/CTBE/splicedeep_AcSS_output_SP803280_D_00-001.gtfacceptor_at_intron.fa.splice2deep  -p /home/bia/sugarcane_introns_local/spliceator_anaysis/BothGenomes_spliceator_acceptor2.fa     -t/home/bia/sugarcane_introns_local/spliceator_anaysis/BothGenomes_spliceator_donor2.fa  -i data/introns/CTBE_introns_IMEter/SP803280_D_All_introns_imeter.fa -b Splice2Deep_analysis/output_s2d/CTBE/splicedeep_DoSS_output_SP803280_D_00-001.gtfdonor_oriza_intron.fa.splice2deep -c Splice2Deep_analysis/output_s2d/CTBE/splicedeep_AcSS_output_SP803280_D_00-001.gtfacceptor_oriza_intron.fa.splice2deep