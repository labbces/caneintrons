# imports
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
