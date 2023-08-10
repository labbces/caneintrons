import re
import csv
file = '/home/bia/sugarcane_introns_local/Souza.final_table.csv'


with open(file, "r") as final_table:
    for line in final_table:
        header = line.split(',')[1]
        if header != "ID":
            k = re.match(r'(.*)\.([0-9]*)-([0-9]*)-([+|-])-(.*)', header)
            full_id = f'{k.group(1)}-{k.group(5)}'
            identif = k.group(1)
            c1 = k.group(2)
            c2 = k.group(3)
            strand = k.group(4)
            event = "INTRON"
            print(full_id, identif, c1, c2, strand, event, sep='\t')
