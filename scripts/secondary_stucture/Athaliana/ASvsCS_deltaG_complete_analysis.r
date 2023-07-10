library(ggplot2)


# How data files were generated:
# 1. Extract SS flanking sequences with extractor 
# 2. Run RNAfold secondary structure prediction 
# 3. Run RNAfold analysis to get seconary stucture features - deltaG.
# - Run for retained data only, constitutive data only and both labeled data


### ANALYSIS OF A. THALIANA DATA ###

# Loading data - retained and constitutive secondary structure databases
#CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_cons_ALL_constitutive_Acceptor_intron.csv')
#CI_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_cons_ALL_constitutive_Donor_intron.csv')
CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/AtRTD2_ALL_constitutive_Acceptor_intron.csv')
CI_donor =read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/AtRTD2_ALL_constitutive_Donor_intron.csv')
CI_acceptor_quasi = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/AtRTDv2_QUASI_ALL_constitutive_Acceptor_intron.csv')
CI_donor_quasi = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/AtRTDv2_QUASI_ALL_constitutive_Donor_intron.csv')


altFive_correct = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_ALTdonor.csv')
altThree_correct = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_ALTacceptor.csv')

retainedIntron_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_INTacceptor.csv')
retainedIntron_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_INTdonor.csv')

cassette_rec_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXrec_donor.csv')
cassette_sk_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXsk_acceptor.csv')
cassette_sk_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXsk_donor.csv')
cassette_rec_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXrec_acceptor.csv')

### Generating histograms ###
# altFive
ggplot(altFive_correct, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 5' Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# altThree
ggplot(altThree_correct, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 5' Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# Intron retention
ggplot(retainedIntron_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Retained intron Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
ggplot(retainedIntron_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Retained intron Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# Cassette 
ggplot(cassette_rec_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Cassette - recognized donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

ggplot(cassette_sk_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Cassette - skipped donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

ggplot(cassette_rec_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Cassette - recognized acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

ggplot(cassette_sk_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Cassette - skipped acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# Constitutive - Donor
ggplot(CI_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") +
  xlim(-100, 5)

ggplot(CI_donor_quasi, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") +
  xlim(-100, 5)
# Constitutive - Acceptor
ggplot(CI_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

ggplot(CI_acceptor_quasi, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")


### Balanced data AS and CI ###
porcentagens = prop.table(table(retainedIntron_donor$deltaG)) * 100
dados_hist = data.frame(delta_g = as.numeric(names(porcentagens)), porcentagens)

porcentagens2 = prop.table(table(CI_donor_quasi$deltaG)) * 100
dados_hist2 = data.frame(delta_g = as.numeric(names(porcentagens2)), porcentagens2)

# ggplot(dados_hist, aes(x = delta_g, fill = "Porcentagem de Ocorrência")) +
#   geom_histogram(position = "identity", alpha = 0.5, color = "#619CFF", fill = '#619CFF', bins = 30) +
#   labs(title = "Stability distributions - Constitutive Acceptor",
#        x = "Valor do Delta G",
#        y = "Porcentagem de Ocorrência") + xlim(-150, 30)


ggplot() +
  geom_bar(data = dados_hist, aes(x = delta_g, y = porcentagens), stat = "identity", alpha = 1, color = "lightblue", fill = "lightblue") +
  geom_bar(data = dados_hist2, aes(x = delta_g, y = porcentagens2), stat = "identity", alpha = 0.1, color = "violet", fill = "violet") +
  labs(title = "Stability distributions - Constitutive Acceptor",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência") +
  xlim(-100, 5)
  


# entre 1%-3% dos arquivos originais...
altFive_correct = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_ALTdonor_balanced.csv')
altThree_correct = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_ALTacceptor_balanced.csv')

retainedIntron_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_INTacceptor_balanced.csv')
retainedIntron_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_INTdonor_balanced.csv')

cassette_rec_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXrec_donor_balanced.csv')
cassette_sk_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXsk_acceptor_balanced.csv')
cassette_sk_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXsk_donor_balanced.csv')

cassette_rec_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_EXrec_acceptor_balanced.csv')

# alternative Five
CI_altFive_constitutive = rbind(CI_acceptor, altFive_correct)
ggplot(CI_altFive_constitutive, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Acceptor vs Alt5 correct)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# alternative Three
CI_altThree_constitutive = rbind(CI_donor, altThree_correct)
ggplot(CI_altThree_constitutive, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Donor x Alt3 correct",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# IR
CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_cons_ALL_constitutive_Acceptor_intron.csv')
CI_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_cons_ALL_constitutive_Donor_intron.csv')
retainedIntron_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_INTacceptor_balanced.csv')
retainedIntron_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_INTdonor_balanced.csv')

CI_IR_donor = rbind(retainedIntron_donor, CI_acceptor)
ggplot(CI_IR_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_IR_acceptor = rbind(retainedIntron_acceptor, CI_donor)
ggplot(CI_IR_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")



CI_IR_donor = rbind(CI_donor, retainedIntron_donor)
ggplot(CI_IR_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

retainedIntron_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/csv/TAIR10_ALL_INTdonor_balanced2.csv')
CI_IR_donor = rbind(CI_donor, retainedIntron_donor)
ggplot(CI_IR_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# Cassette
CI_cassette_donor = rbind(CI_donor, cassette_rec_donor)
ggplot(CI_cassette_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exons - Recognized donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_cassette_donor = rbind(CI_donor, cassette_sk_donor)
ggplot(CI_cassette_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exon - Skipped exon donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_cassette_acceptor = rbind(CI_acceptor, cassette_rec_acceptor)
ggplot(CI_cassette_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exons - Recognized acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_cassette_acceptor = rbind(CI_acceptor, cassette_sk_acceptor)
ggplot(CI_cassette_aceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exon - Skipped acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")


### T-TEST ###

# Constitutive splicing vs alternative 3'
t.test(CI_donor$deltaG, altThree_correct$deltaG) 
# Constitutive splicing vs alternative 5'
t.test(CI_acceptor$deltaG, altFive_correct$deltaG) 
# Constitutive splicing vs retined intron (Donor)
t.test(CI_donor$deltaG, retainedIntron_donor$deltaG) 
# Constitutive splicing vs retined intron (Acceptor)
t.test(CI_acceptor$deltaG, retainedIntron_acceptor$deltaG) 

# Constitutive splicing vs recognized acceptor 
t.test(CI_acceptor$deltaG, cassette_rec_acceptor$deltaG) 
# Constitutive splicing vs  skipped acceptor
t.test(CI_acceptor$deltaG, cassette_sk_acceptor$deltaG) 
# Recognized Acceptor vs skipped acceptor
t.test(cassette_rec_acceptor$deltaG, cassette_sk_acceptor$deltaG) 

# Constitutive splicing vs recognized donor 
t.test(CI_donor$deltaG, cassette_rec_donor$deltaG) 
# Constitutive splicing vs skipped donor
t.test(CI_donor$deltaG, cassette_sk_donor$deltaG) 
# Recognized Donor vs skipped donor
t.test(cassette_rec_donor$deltaG, cassette_sk_donor$deltaG) 



