import pandas as pd
#file = '/home/bia/sugarcane_introns_local/NewIntrons_STAR/CTBE_NewIntrons_STAR.csv'
#output = '/home/bia/sugarcane_introns_local/NewIntrons_STAR/CTBE_NewIntrons_STAR_filtered.csv'

file = '/home/bia/sugarcane_introns_local/NewIntrons_STAR/Souza_NewIntrons_STAR.csv'
output = '/home/bia/sugarcane_introns_local/NewIntrons_STAR/Souza_NewIntrons_STAR_filtered.csv'

df = pd.read_csv(file)
q1 = df['mapped_sum'].quantile(0.25)
print("Q1 value:", q1)
filtered_df = df[df['mapped_sum'] >= q1]
filtered_df.to_csv(output, index=False)
