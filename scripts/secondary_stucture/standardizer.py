import re
import pandas as pd
### HUMAN ###
# Constitutive coordinates - output of selectConstitutiveExons.py
'''
cs = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/Human_genome/generated/Cexons_coord_human_hg18.txt', sep='\t', header=None)
cs['filt'] = cs[[1, 2, 3, 4]].astype(str).apply('_'.join, axis=1)
cs = cs.drop_duplicates('filt')
cs = cs.drop('filt', axis=1)
print('human')
print(cs)
cs.to_csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated/Cexons_coord_human_hg18.tsv',
          sep='\t', index=False, header=None)
'''

# Alternative splicing - from KnownAlt
'''
as_f = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/Human_genome/downloaded/knownAlt_hg18.txt', sep="\t", header=None)
as_f = as_f[[0, 1, 2, 3, 6, 4]]
as_f = as_f.replace({
    'retainedIntron': 'INT',
    'cassetteExon': 'EX',
    'altFivePrime': 'ALTA',
    'altThreePrime': 'ALTD'
})
as_f.to_csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated/knownAlt_hg18.tsv',
            sep="\t", header=None, index=False)
'''

### ARABIDOPSIS ###
# Constitutive coordinates - output of selectConstitutiveExons.py
'''
cs = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/Cexons_coord_athaliana_AtRTDv2_QUASI.tsv', sep='\t', header=None)
cs['filt'] = cs[[1, 2, 3, 4]].astype(str).apply('_'.join, axis=1)
cs = cs.drop_duplicates('filt')
cs = cs.drop('filt', axis=1)
print('arabidopsis')
cs.to_csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/Cexons_coord_athaliana_AtRTDv2_QUASI.tsv',
          sep='\t', index=False, header=None)
'''

# Alternative splicing - from EVENT
# cut -f 2,3,7 EVENT_INFO-araTha10.tab >> EVENT_INFO-araTha10_simplified.tab
'''
as_file = '/home/bia/sugarcane_introns_local/data/Athaliana_genome/downloaded/EVENT_INFO-araTha10_simplified.tab'
with open(as_file, 'r') as as_file:
    for line in as_file:
        id, coord, strand = line.split('\t')
        chr, coord = coord.split(':')
        c1, c2 = coord.split('-')
        strand = strand.split(':')[2].strip()
        e = re.match(r'Ath([A-Z]*)[0-9]*', id)
        event = e.group(1).strip()
        # print(event)
        print(id, chr, c1, c2, strand, event, sep='\t')
as_f = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/EVENT_INFO-araTha10_simplified.txt', sep='\t', header=None)
as_f['filt'] = as_f[[1, 2, 3, 4]].astype(str).apply('_'.join, axis=1)
as_f = as_f.drop_duplicates('filt')
as_f = as_f.drop('filt', axis=1)
print(as_f)
as_f.to_csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/EVENT_INFO-araTha10_simplified.tsv',
            sep='\t', index=False, header=None)
'''


### MAIZE ###
# Constitutive coordinates - output of selectConstitutiveExons.py
'''
cs = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/maize/generated/Cexons_coord_maize_ysu2017.txt', sep='\t', header=None)
cs['filt'] = cs[[1, 2, 3, 4]].astype(str).apply('_'.join, axis=1)
cs = cs.drop_duplicates('filt')
cs = cs.drop('filt', axis=1)
print('maize')
print(cs)
cs.to_csv('/home/bia/sugarcane_introns_local/data/maize/generated/Cexons_coord_maize_ysu2017.tsv',
          sep='\t', index=False, header=None)
'''
# Alternative splicing - from landscape
as_f = pd.read_csv(
    '/home/bia/sugarcane_introns_local/data/maize/downloaded/landscape.tab', sep="\t", header=None)
as_f = as_f[[1, 0, 3, 4, 5, 7]]

as_f = as_f.replace({
    'intron retention': 'INT',
    'exon skipping': 'EX',
    'alternative acceptor site': 'ALTA',
    'alternative donor site': 'ALTD'
})
print(as_f)
as_f.to_csv('/home/bia/sugarcane_introns_local/data/maize/downloaded/landscape.tsv',
            sep="\t", header=None, index=False)
