library(ggplot2)
library(randomForest)

# How data files were generated:
# 1. Extract SS flanking sequences with extractor 
# 2. Run RNAfold secondary structure prediction 
# 3. Run RNAfold analysis to get seconary stucture features - deltaG.
# - Run for retained acceptor and donor, constitutive data for acceptor and donor

### Loading data ###
CI_donor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Donor1.csv')
RI_donor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/retained_Athaliana_Donor1.csv')
CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Acceptor1.csv')
RI_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/retained_Athaliana_Acceptor1.csv')
data = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/general_data_1.csv')

### Violin plots to visualize data dispersion ###
ggplot(data, aes(x=label, y=freq, fill=label)) + geom_violin() +  labs(title="Paired nucleotides frequency")
ggplot(data, aes(x=label, y=freq1.2, fill=label)) + geom_violin() + labs(title="Paired nucleotides frequency - 1/2 sequence")
ggplot(data, aes(x=label, y=freq2.2, fill=label)) + geom_violin()  + labs(title="Paired nucleotides frequency - 2/2 sequence")
ggplot(data, aes(x=label, y=freq1.3, fill=label)) + geom_violin() + labs(title="Paired nucleotides frequency - 1/3 sequence")
ggplot(data, aes(x=label, y=freq2.3, fill=label)) + geom_violin() + labs(title="Paired nucleotides frequency - 2/3 sequence")
ggplot(data, aes(x=label, y=freq3.3, fill=label)) + geom_violin() + labs(title="Paired nucleotides frequency - 3/3 sequence")
ggplot(data, aes(x=label, y=substru, fill=label)) + geom_violin()  + labs(title="Substrutures Amount")
ggplot(data, aes(x=label, y=deltaG, fill=label)) + geom_violin() + labs(title="Full sequence Free energy")
ggplot(data, aes(x=label, y=Stem_comp, fill=label)) + geom_violin() + labs(title="Full sequence Average stem length")


### T tests to identify relevant differences in data ###
# Paired nucleotides frequency for donor and acceptor
t.test(CI_donor$freq, RI_donor$freq) 
t.test(CI_acceptor$freq, RI_acceptor$freq) 
# Paired nucleotides frequency for donor and acceptor 1/2
t.test(CI_donor$freq1.2, RI_donor$freq1.2) 
t.test(CI_acceptor$freq1.2, RI_acceptor$freq1.2) 
# Paired nucleotides frequency for donor and acceptor 2/2
t.test(CI_donor$freq2.2, RI_donor$freq2.2) 
t.test(CI_acceptor$freq2.2, RI_acceptor$freq2.2) 
# Paired nucleotides frequency for donor and acceptor 1/3
t.test(CI_donor$freq1.3, RI_donor$freq1.3) 
t.test(CI_acceptor$freq1.3, RI_acceptor$freq1.3) 
# Paired nucleotides frequency for donor and acceptor 2/3
t.test(CI_donor$freq2.3, RI_donor$freq2.3) 
t.test(CI_acceptor$freq2.3, RI_acceptor$freq2.3) 
# Paired nucleotides frequency for donor and acceptor 3/3
t.test(CI_donor$freq3.3, RI_donor$freq3.3) 
t.test(CI_acceptor$freq3.3, RI_acceptor$freq3.3) 
# Substrutures Amount for donor and acceptor
t.test(CI_donor$substru, RI_donor$substru) 
t.test(CI_acceptor$substru, RI_acceptor$substru) 
# Free energy for donor and acceptor
t.test(CI_donor$deltaG, RI_donor$deltaG) 
t.test(CI_acceptor$deltaG, RI_acceptor$deltaG) 
# Stem length for donor and acceptor
t.test(CI_donor$Stem_comp, RI_donor$Stem_comp) 
t.test(CI_acceptor$Stem_comp, RI_acceptor$Stem_comp) 

### RanfomForest ML###
balanced_data = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/balanced_Athaliana_certo.csv')
balanced_data$label = factor(balanced_data$label)
randomForest(label ~ .,data=balanced_data, ntree=500)


