library ('data.table')
library(ggplot2)
# How data files were generated:
# 1. Extract SS flanking sequences with extractor 
# 2. Run RNAfold secondary structure prediction 
# 3. Run RNAfold analysis to get seconary stucture features - deltaG.
# - Run for retained data only, constitutive data only and both labeled data

# Loading data - retained and constitutive secondary structure databases
CI_donor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Donor1.csv')
AS_donor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/as_Athaliana_Donor.csv')
CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Acceptor1.csv')
AS_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/as_Athaliana_Acceptor.csv')


### Generating histograms ###
# AS - Donor
ggplot(AS_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="#619CFF",color="#619CFF") +
  labs(title="Stability distributions - Alternative splicing Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") 

# AS - Acceptor
ggplot(AS_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="#619CFF",color="#619CFF") +
  labs(title="Stability distributions - Alternative splicing Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") 


# Constitutive - Donor
ggplot(CI_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="pink") +
  labs(title="Stability distributions - Constiitutive Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# Constitutive - Acceptor
ggplot(CI_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="pink") +
  labs(title="Stability distributions - Constiitutive Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

### Balanced data RI and CI ###
data1 = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/asci_Athaliana_donor.csv')
ggplot(data1, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

data2 = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/asci_Athaliana_acceptor.csv')
ggplot(data2, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")
