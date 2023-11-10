import pandas as pd
file = "/home/bia/sugarcane_introns_local/data/RNAseq_output/CTBE_SJ_comSRR.out.tab"
#file = "/home/bia/sugarcane_introns_local/data/RNAseq_output/Souza_SJ_comSRR.out.tab"

#output = "/home/bia/sugarcane_introns_local/data/newIntrons/CTBE_NewIntrons_STAR"
#output = "/home/bia/sugarcane_introns_local/data/newIntrons/Souza_NewIntrons_STAR"
output = "xxxxxx"

# STAR output dict
full_data = {}
noted_data = {}
# Chromosome - Start - Enr - Strand - [Unmap, MultiMap, Sum]

# Abre o STAR output
with open(file, 'r') as files:
    for line in files:
        chr, start, end, strand, canonical, noted, un_map, multi_map, overhang, srr = line.split(
            '\t')
        start = int(start)
        end = int(end)
        noted = int(noted)
        un_map = int(un_map)
        multi_map = int(multi_map)
        strand = int(strand)
        # print(chr, start, end, strand, canonical,
        # noted, un_map, multi_map, overhang, srr)
        # Filtra pelo tamanho
        if strand != 0:  # remove undefined strands
            # apenas introns com tamanho apropriado
            # if ((int(end) - int(start)) + 1) >= 100 and ((int(end) - int(start)) + 1) <= 150:
            # filtra para ficar apenas com introns não anotados
            if noted == 1:
                if chr not in noted_data.keys():
                    noted_data[chr] = {}
                    noted_data[chr][start] = {}
                    noted_data[chr][start][end] = 1
                else:
                    if start not in noted_data[chr].keys():
                        noted_data[chr][start] = {}
                        noted_data[chr][start][end] = 1
                    else:
                        if end not in noted_data[chr][start].keys():
                            noted_data[chr][start][end] = 1
with open(file, 'r') as files:
    for line in files:
        chr, start, end, strand, canonical, noted, un_map, multi_map, overhang, srr = line.split(
            '\t')
        start = int(start)
        end = int(end)
        noted = int(noted)
        un_map = int(un_map)
        multi_map = int(multi_map)
        strand = int(strand)
        # print(chr, start, end, strand, canonical,
        # noted, un_map, multi_map, overhang, srr)
        # Filtra pelo tamanho
        if strand != 0:  # remove undefined strands
            # if ((int(end) - int(start)) + 1) >= 100 and ((int(end) - int(start)) + 1) <= 150:
            if noted == 0:  # 0 = no annotation 1: with annotation
                # print(full_data)
                if chr not in noted_data.keys():
                    if chr not in full_data.keys():  # novo cromossomo - add chr, start, end, info
                        full_data[chr] = {}
                        full_data[chr][start] = {}
                        full_data[chr][start][end] = {}
                        full_data[chr][start][end][strand] = [
                            un_map, multi_map, (int(un_map)+int(multi_map))]
                    else:
                        test_start = False  # para avaliar se start ja foi add em full_data
                        # criando intervalo de start - start+5 start-5
                        for s in range(start-5, start+6):
                            # check se start já foi add
                            if s in full_data[chr].keys():
                                test_start = True  # se start já tiver sido add, test_start = True

                                test_end = False  # para avaliar se end ja foi add em full_data
                                # criando intervalo de end - end+5 end-5
                                for e in range(end-5, end+6):
                                    if e in full_data[chr][s].keys():
                                        test_end = True  # se end já tiver sido add, test_end = True
                                        # add info
                                        if strand not in full_data[chr][s][e].keys():
                                            full_data[chr][s][e][strand] = [
                                                un_map, multi_map, (int(un_map)+int(multi_map))]
                                        else:
                                            full_data[chr][s][e][strand][0] += un_map
                                            full_data[chr][s][e][strand][1] += multi_map
                                            full_data[chr][s][e][strand][2] += (
                                                un_map+multi_map)
                                        break
                                if test_end == False:  # end não foi adicionado - test_end = False
                                    full_data[chr][s][end] = {}
                                    full_data[chr][s][end][strand] = [
                                        un_map, multi_map, (int(un_map)+int(multi_map))]
                                break
                        if test_start == False:  # start não foi adicioando - test_start = False
                            # add start, end e info
                            full_data[chr][start] = {}
                            full_data[chr][start][end] = {}
                            full_data[chr][start][end][strand] = [un_map, +
                                                                  multi_map, (int(un_map)+int(multi_map))]
                else:
                    test = False
                    for s in range(start-5, start+6):
                        if s in noted_data[chr].keys():
                            if any(e in noted_data[chr][s].keys() for e in range(end-5, end+6)):
                                test = True
                                break
                    if test == False:
                        for s in range(start-5, start+6):
                            test_start = False
                            if chr in full_data.keys():
                                # check se start já foi add
                                if s in full_data[chr].keys():
                                    test_start = True  # se start já tiver sido add, test_start = True

                                    test_end = False  # para avaliar se end ja foi add em full_data
                                    # criando intervalo de end - end+5 end-5
                                    for e in range(end-5, end+6):
                                        if e in full_data[chr][s].keys():
                                            test_end = True  # se end já tiver sido add, test_end = True
                                            # add info
                                            if strand not in full_data[chr][s][e].keys():
                                                full_data[chr][s][e][strand] = [
                                                    un_map, multi_map, (int(un_map)+int(multi_map))]
                                            else:
                                                full_data[chr][s][e][strand][0] += un_map
                                                full_data[chr][s][e][strand][1] += multi_map
                                                full_data[chr][s][e][strand][2] += (
                                                    un_map+multi_map)
                                            break
                                    if test_end == False:  # end não foi adicionado - test_end = False
                                        full_data[chr][s][end] = {}
                                        full_data[chr][s][end][strand] = [
                                            un_map, multi_map, (int(un_map)+int(multi_map))]
                                    break
                                if test_start == False:  # start não foi adicioando - test_start = False
                                    # add start, end e info
                                    full_data[chr][start] = {}
                                    full_data[chr][start][end] = {}
                                    full_data[chr][start][end][strand] = [un_map, +
                                                                          multi_map, (int(un_map)+int(multi_map))]
                            else:
                                full_data[chr] = {}
                                full_data[chr][start] = {}
                                full_data[chr][start][end] = {}
                                full_data[chr][start][end][strand] = [
                                    un_map, multi_map, (int(un_map)+int(multi_map))]


'''total_introns = 0
for chr in full_data.keys():
    # print(len(full_data))  # quantidade de chromosomos/leituras
    # print(len(full_data[chr]))  # quantidade de starts
    for start in full_data[chr].keys():
        for end in full_data[chr][start].keys():
            for strand in full_data[chr][start][end].keys():
                if full_data[chr][start][end][strand][2] >= 62:
                    total_introns += 1
print(total_introns)'''


strand_db = {1: "+", 2: "-"}
with open(f'{output}_nonFiltered.csv', 'w') as out:
    out.write('contig,start,end,strand,uniquely_mapped,multi_mapped,mapped_sum\n')
    for id in full_data.keys():
        for start in full_data[id].keys():
            for end in full_data[id][start].keys():
                for strand, data in full_data[id][start][end].items():
                    strand = strand_db[strand]
                    datas = list(map(str, data))
                    datass = [str(item) for item in datas]
                    datas = ','.join(datass)
                    out.write(f'{id},{start},{end},{strand},{datas}\n')

df = pd.read_csv(f'{output}_nonFiltered.csv')
#df = df[df['mapped_sum'] >= 62]
df.to_csv(f"{output}_filtered.csv", index=False)
