
# How data files were generated:
# 1. Extract SS flanking sequences with extractor 
# 2. Run RNAfold secondary structure prediction 
# 3. Run RNAfold analysis to get seconary stucture features - deltaG.
# - Run for retained data only, constitutive data only and both labeled data

# Loading data - retained and constitutive secondary structure databases
data_table = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/retained_balanced_Athaliana.csv')
data_table2 = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana.csv')


### Generating histograms ###
# Retained - Donor
ggplot(data_table, aes(x=deltaG_d, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="#619CFF",color="#619CFF") +
  labs(title="Stability distributions - Reteined Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") 

# Retained - Acceptor
ggplot(data_table, aes(x=deltaG_a, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="#619CFF",color="#619CFF") +
  labs(title="Stability distributions - Reteined Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") 


# Constitutive - Donor
ggplot(data_table2, aes(x=deltaG_d, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="pink") +
  labs(title="Stability distributions - Constiitutive Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# Constitutive - Acceptor
ggplot(data_table2, aes(x=deltaG_a, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="pink") +
  labs(title="Stability distributions - Constiitutive Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

### Balanced data RI and CI ###
df = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/full_balanced_Athaliana.csv')
ggplot(df, aes(x=deltaG_a, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")



  