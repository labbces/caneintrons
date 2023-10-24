import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2_unweighted


headers = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/newIntrons/fastas/Souza_newIntrons_Splice2Deep_acceptor_intron.fa', header=None)
donor_at = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/newIntrons/s2d_output/splicedeep_DoSS_output_Souza_at.splice2deep', header=None)
acceptor_at = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/newIntrons/s2d_output/splicedeep_AcSS_output_Souza_at.splice2deep', header=None)
donor_or = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/newIntrons/s2d_output/splicedeep_DoSS_output_Souza_oriza.splice2deep', header=None)
acceptor_or = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/newIntrons/s2d_output/splicedeep_AcSS_output_Souza_oriza.splice2deep', header=None)

headers = headers[headers.iloc[:, 0].str.startswith('>')]
headers = headers.reset_index(drop=True)

df = pd.concat([headers, donor_at, acceptor_at,
                donor_or, acceptor_or], axis=1)
df.columns = ['ID', 'Donor_at', 'Acceptor_at',
              'Donor_Oriza', 'Acceptor_Oriza']

print(df)
#df.to_csv('Souza_newIntrons_s2dAnalysis.csv', index=False, header=False, sep='\t')

good_good = df[(df['Donor_at'] == 1) & (df['Acceptor_at'] == 1) & (
    df['Donor_Oriza'] == 1) & (df['Acceptor_Oriza'] == 1)]

# good_introns: Linhas onde Donor_at e Acceptor_at sejam ambos igual a 1 OU onde Donor_Oriza e Acceptor_Oriza sejam iguais a 1
good_introns = df[((df['Donor_at'] == 1) & (df['Acceptor_at'] == 1)) | (
    (df['Donor_Oriza'] == 1) & (df['Acceptor_Oriza'] == 1))]

# good_at: Linhas onde Donor_at e Acceptor_at sejam ambos igual a 1
good_at = df[(df['Donor_at'] == 1) & (df['Acceptor_at'] == 1)]

# good_oriza: Linhas onde Donor_Oriza e Acceptor_Oriza sejam ambos igual a 1
good_oriza = df[(df['Donor_Oriza'] == 1) & (df['Acceptor_Oriza'] == 1)]

print(good_good)
print(good_introns)
print(good_at)
print(good_oriza)

ga = good_at.shape[0] - good_good.shape[0]
go = good_oriza.shape[0] - good_good.shape[0]
int = good_good.shape[0]
good_introns.to_csv('Souza_newIntrons_s2dAnalysis.csv',
                    index=False, header=True, sep='\t')
v = venn2_unweighted(subsets=(go, ga, int),
                     set_labels=('Oriza', 'Arabidopsis'))
v.get_patch_by_id('10').set_color('blue')  # Conjunto A
v.get_patch_by_id('01').set_color('pink')  # Conjunto B
v.get_patch_by_id('11').set_color('purple')
# Adicione um título
plt.title("Souza")

# Exiba o gráfico
plt.show()

'''# imports
import glob
import matplotlib.pyplot as plt
from matplotlib_venn import venn2_unweighted, venn2_circles

all_data = {'shared_ids_at': 0, 'good_introns_at': 0, 'shared_ids_Or': 0, 'good_introns_Or': 0, 'shared_ids_Both': 0, 'good_introns_Both': 0}

  
compfile_list = glob.glob('*_Souza-compFile.txt')

for filename in compfile_list:
    with open(filename) as f:
        lines = f.readlines()

    shared_ids_at = int(lines[0].split(': ')[1])
    all_data['shared_ids_at'] = all_data['shared_ids_at'] + shared_ids_at
    good_introns_at = int(lines[1].split(': ')[1])
    all_data['good_introns_at'] = good_introns_at + all_data['good_introns_at']

    shared_ids_or = int(lines[2].split(': ')[1])
    all_data['shared_ids_Or'] = all_data['shared_ids_Or'] + shared_ids_or
    good_introns_or = int(lines[3].split(': ')[1])
    all_data['good_introns_Or'] = good_introns_or + all_data['good_introns_Or']

    shared_ids = int(lines[4].split(': ')[1])
    all_data['shared_ids_Both'] = all_data['shared_ids_Both'] + shared_ids
    good_introns = int(lines[5].split(': ')[1])
    all_data['good_introns_Both'] = good_introns + all_data['good_introns_Both']


with open('Souza_splice2deep_analysis.txt', 'w') as f:
#with open('CTBE_splice2deep_analysis.txt', 'w') as f:
    f.write(f'Number of shared ids - spliceator and splice2deep at model: {all_data["shared_ids_at"]}\n')
    f.write(f'Number of good introns - splice2deep at model : {all_data["good_introns_at"]}\n')
    f.write(f'Number of shared ids - spliceator and splice2deep oriza model: {all_data["shared_ids_Or"]}\n')
    f.write(f'Number of good introns - splice2deep oriza model: {all_data["good_introns_Or"]}\n')
    f.write(f'Number of shared ids - spliceator and splice2deep both models: {all_data["shared_ids_Both"]}\n')
    f.write(f'Number of good introns - splice2deep both models: {all_data["good_introns_Both"]}\n')

def create_venn_diagram(data, name, label_group):
    plt.figure()
    venn = venn2_unweighted(subsets=data, set_labels=label_group)

    venn.get_label_by_id('10').set_text(f'{data[0]}')
    venn.get_label_by_id('01').set_text(f'{data[1]}')
    venn.get_label_by_id('11').set_text(f'{data[2]}')

    plt.savefig(name)

data_good_introns = [ all_data["good_introns_at"] - all_data["good_introns_Both"], all_data["good_introns_Or"] - all_data["good_introns_Both"] ,all_data["good_introns_Both"]]
create_venn_diagram(data_good_introns, 'Souza_good_introns_VennD.png', ('Arabidopsis', 'Oriza'))

good_introns_splicator = 36277
# ctbe = 11026
# IQ USP = 36277
# data_shaeredIDs = [good_introns_splicator-all_data["shared_ids_Both"], all_data["good_introns_Both"] - all_data["shared_ids_Both"], all_data["shared_ids_Both"]]
data_shaeredIDs = [good_introns_splicator-all_data["shared_ids_at"], all_data["good_introns_at"] - all_data["shared_ids_at"], all_data["shared_ids_at"]]
create_venn_diagram(data_shaeredIDs, 'Souza_shared_introns_at_VennD.png', ('Splcieator', 'Splice2deep'))
'''
