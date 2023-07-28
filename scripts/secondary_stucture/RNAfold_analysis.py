#TODO if run extractor again, get event auto by header pattern
import pandas as pd
import re
import argparse
import subprocess


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True,
                    help='RNAfold output file')
parser.add_argument("-o", "--output", type=str, required=True,
                    help='For "teste" as -o, the output will be teste.csv')
#parser.add_argument("-t", "--tag", type=str, required=True)
parser.add_argument("-e", "--event", type=str, required=True,
                    help='ALTA, INT, ALTD, CS or EX')
parser.add_argument("-g", "--grapple", type=str, required=False,
                    help='path to GraPPLE get_properties.R script - if avaiable')
args = parser.parse_args()


# Calling grapple get properties R script as a function
def grapple_prop(script, structure):
    output = subprocess.check_output(['Rscript', script, structure])
    output = output.decode('utf-8')
    return output


# infile = '/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Acceptor1.RNAfold'
# outfile = '/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Acceptor1.csv'
infile = args.input
outfile = args.output

with open(infile, 'r') as infile, open(f"{outfile}.csv", 'w') as outfile:
    if args.grapple is not None:
        outfile.write(
            f'id,label,freq,freq1/2,freq2/2,freq1/3,freq2/3,freq3/3,substru,deltaG,Stem_comp,articulation_points,average_path_length,average_node_betwennness,variance_node_betweenness,average_edge_betweenness,variance_edge_betweenness,average_cocitation_coupling,average_bibliographic_coupling,average_closeness_centrality,variance_closeness_centrality,average_burts_constraint,variance_burts_constraint,average_degree,diameter,girth,average_coreness,varince_coreness,max_coreness,graph_density,transitivity\n')
    else:
        outfile.write(
            'id,label,freq,freq1/2,freq2/2,freq1/3,freq2/3,freq3/3,substru,deltaG,Stem_comp\n')

    for line in infile:
        line = line.rstrip('\n')
        # print(line)
        if line.startswith('>'):
            id = line
        elif line.rstrip().endswith(')'):
            struc = line.split(' ', 1)[0]
            deltaG = float(line.split(' ', 1)[1].lstrip('(').rstrip(')'))
            label = args.event

            # calculando a frequencia geral
            comprim = len(struc)
            freq = (struc.count('(') + struc.count(')'))/comprim

            # calculando posição 1/3
            comprim3 = int(comprim/3)

            struc13 = struc[0:comprim3]
            struc23 = struc[comprim3:(comprim3*2)]
            struc33 = struc[(comprim3*2):]

            freq13 = (struc13.count('(') + struc13.count(')'))/comprim3
            freq23 = (struc23.count('(') + struc23.count(')'))/comprim3
            freq33 = (struc33.count('(') + struc33.count(')'))/comprim3

            # calculando posição 1/2
            comprim2 = int(comprim/2)

            struc12 = struc[0:comprim2]
            struc22 = struc[comprim2:]

            freq12 = (struc12.count('(') + struc12.count(')'))/comprim2
            freq22 = (struc22.count('(') + struc22.count(')'))/comprim2

            # calculando sub estruturas
            substru = re.findall(r'\)\.*\(', struc)
            substru = len(substru) + 1
            if freq == 0:
                substru = 0

            # calculando comprimento do stem
            n = re.findall(r'(\(\(*\()', struc)
            o = re.findall(r'(\)\)*\))', struc)
            comprim = len(n) + len(o)
            semistem_comp = 0
            if comprim > 0:
                for stem in n:
                    semistem_comp += len(stem)
                for stem in o:
                    semistem_comp += len(stem)
                stem_comp = semistem_comp/comprim
            else:
                stem_comp = 0

            if args.grapple is not None:
                out_grapple = grapple_prop(
                    script=args.grapple, structure=struc)
                outfile.write(
                    f'{id},{label},{freq},{freq12},{freq22},{freq13},{freq23},{freq33},{substru},{deltaG},{stem_comp},{out_grapple.strip()}\n')
            else:
                outfile.write(
                    f'{id},{label},{freq},{freq12},{freq22},{freq13},{freq23},{freq33},{substru},{deltaG},{stem_comp}\n')
'''

df1 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/retained_Athaliana_Donor1.csv')
columns = ['label_d', 'freq_d', 'freq1/2_d', 'freq2/2_d', 'freq1/3_d',
           'freq2/3_d', 'freq3/3_d', 'substru_d', 'deltaG_d', 'Stem_comp_d']
df1.columns = columns

df2 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/retained_Athaliana_Acceptor1.csv')
columns = ['label_a', 'freq_a', 'freq1/2_a', 'freq2/2_a', 'freq1/3_a',
           'freq2/3_a', 'freq3/3_a', 'substru_a', 'deltaG_a', 'Stem_comp_a']
df2.columns = columns
# print(df1)
# print(df2)

#df = df1.merge(df2, on=index, how='inner')
df = df1.join(df2)
df = df.drop(columns=['label_d', 'label_a'])
df['label'] = 'retained'
# print(df)

df1 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/balanced_constitutive_donor1.csv')
columns = ['label_d', 'freq_d', 'freq1/2_d', 'freq2/2_d', 'freq1/3_d',
           'freq2/3_d', 'freq3/3_d', 'substru_d', 'deltaG_d', 'Stem_comp_d']
df1.columns = columns
df2 = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/balanced_constitutive_acceptor1.csv')
columns = ['label_a', 'freq_a', 'freq1/2_a', 'freq2/2_a', 'freq1/3_a',
           'freq2/3_a', 'freq3/3_a', 'substru_a', 'deltaG_a', 'Stem_comp_a']
df2.columns = columns
print(df1)
print(df2)

#df = df1.merge(df2, on=index, how='inner')
dff = df1.join(df2)
print(dff)
dff = dff.drop(columns=['label_d', 'label_a'])
dff['label'] = 'constitutive'
print(dff)

df = df.reset_index(drop=True)
df_final = pd.concat([df, dff])
print(df_final)

df_final.to_csv(f'full_balanced_Athaliana.csv', index=False)
df.to_csv(f'retained_balanced_Athaliana.csv', index=False)
dff.to_csv(f'constitutive_balanced_Athaliana.csv', index=False)
'''
