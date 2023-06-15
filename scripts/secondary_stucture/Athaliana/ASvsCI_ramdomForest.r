library(ggplot2)
library(randomForest)

# How data files were generated:
# 1. Extract SS flanking sequences with extractor 
# 2. Run RNAfold secondary structure prediction 
# 3. Run RNAfold analysis to get seconary stucture features - deltaG.
# - Run for retained acceptor and donor, constitutive data for acceptor and donor

### Loading data ###
CI_donor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Donor1.csv')
AS_donor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/as_Athaliana_Donor.csv')
CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/constitutive_Athaliana_Acceptor1.csv')
AS_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/as_Athaliana_Acceptor.csv')
data = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/as_Athaliana_general.csv')

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
t.test(CI_donor$freq, AS_donor$freq) 
t.test(CI_acceptor$freq, AS_acceptor$freq) 
# Paired nucleotides frequency for donor and acceptor 1/2
t.test(CI_donor$freq1.2, AS_donor$freq1.2) 
t.test(CI_acceptor$freq1.2, AS_acceptor$freq1.2) 
# Paired nucleotides frequency for donor and acceptor 2/2
t.test(CI_donor$freq2.2, AS_donor$freq2.2) 
t.test(CI_acceptor$freq2.2, AS_acceptor$freq2.2) 
# Paired nucleotides frequency for donor and acceptor 1/3
t.test(CI_donor$freq1.3, AS_donor$freq1.3) 
t.test(CI_acceptor$freq1.3, AS_acceptor$freq1.3) 
# Paired nucleotides frequency for donor and acceptor 2/3
t.test(CI_donor$freq2.3, AS_donor$freq2.3) 
t.test(CI_acceptor$freq2.3, AS_acceptor$freq2.3) 
# Paired nucleotides frequency for donor and acceptor 3/3
t.test(CI_donor$freq3.3, AS_donor$freq3.3) 
t.test(CI_acceptor$freq3.3, AS_acceptor$freq3.3) 
# Substrutures Amount for donor and acceptor
t.test(CI_donor$substru, AS_donor$substru) 
t.test(CI_acceptor$substru, AS_acceptor$substru) 
# Free energy for donor and acceptor
t.test(CI_donor$deltaG, AS_donor$deltaG) 
t.test(CI_acceptor$deltaG, AS_acceptor$deltaG) 
# Stem length for donor and acceptor
t.test(CI_donor$Stem_comp, AS_donor$Stem_comp) 
t.test(CI_acceptor$Stem_comp, AS_acceptor$Stem_comp) 

### RanfomForest ML###
balanced_data = read.csv('/home/bia/sugarcane_introns_local/data/RNAstruc_Athaliana/pergunta1/as_Athaliana_asci_general.csv')
balanced_data$label = factor(balanced_data$label)
randomForest(label ~ .,data=balanced_data, ntree=500)
help(randomForest)

