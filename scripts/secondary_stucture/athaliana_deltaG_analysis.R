library(ggplot2)

# How data files were generated:
# 1.Get coordinated of AS and CS events (AS and selectConstitutiveExons.py)
# 2.Extract SS flanking sequences with extractor_splicing.py using genome and coordinate file
# 3.Run RNAfold secondary structure prediction 
# 4.Run RNAfold_analysis.py to get seconary stucture features - deltaG.

##########################################
### ANALYSIS FOR Athaliana_genome DATA ###
##########################################

### LOADING CSV FILES - from step 4 ###
# constitutive events
CS_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_CS_donor.csv')
CS_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_CS_acceptor.csv')

# ALTA events
ALTD_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_ALTD_donor.csv')

# ALTD events
ALTA_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_ALTA_acceptor.csv')

# INT events
INT_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_INT_donor.csv')
INT_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_INT_acceptor.csv')

# EX events
EX_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_EXsk_donor.csv')
EX_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_EXsk_acceptor.csv')

### INDIVIDUAL ANALYSIS ### 
# constitutive events
ggplot(CS_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="royalblue",color="royalblue") +
  labs(title="Constitutive splicing - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

ggplot(CS_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="royalblue",color="royalblue") +
  labs(title="Constitutive splicing - Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# ALTA events
ggplot(ALTA_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Alternative Acceptor/5' - Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# ALTD events
ggplot(ALTD_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Alternative Donor/3' - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# INT events
ggplot(INT_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Intron retained - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

ggplot(INT_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Intron retained - Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

# EX events
ggplot(EX_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Exon skipping - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)

ggplot(EX_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Exon skipping - Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)




### COMPARATIVE ANALYSIS ###
#------# DONOR SITE #------#
CS_donor_porcent = prop.table(table(CS_donor$deltaG)) * 100
dados_hist_CS_donor = data.frame(delta_g = as.numeric(names(CS_donor_porcent)), CS_donor_porcent)

# ALTD
ALTD_porcent = prop.table(table(ALTD_donor$deltaG)) * 100
dados_hist_ALTD = data.frame(delta_g = as.numeric(names(ALTD_porcent)), ALTD_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Constitutive Donor"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Stability distributions - Constitutive Donor vs ALTD",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Constitutive Donor" = "royalblue", "ALTD" = "pink"))

# INT donor
INT_donor_porcent = prop.table(table(INT_donor$deltaG)) * 100
dados_hist_INT_donor = data.frame(delta_g = as.numeric(names(INT_donor_porcent)), INT_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Constitutive Donor"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_donor, aes(x = delta_g, y = INT_donor_porcent, fill = "INT donor"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Stability distributions - Constitutive Donor vs INT donor",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Constitutive Donor" = "royalblue", "INT donor" = "pink"))

# EX donor
EX_donor_porcent = prop.table(table(EX_donor$deltaG)) * 100
dados_hist_EX_donor = data.frame(delta_g = as.numeric(names(EX_donor_porcent)), EX_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Constitutive Donor"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX donor"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Stability distributions - Constitutive Donor vs EX donor",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Constitutive Donor" = "royalblue", "EX donor" = "pink"))


#------# ACCEPTOR SITE #------#
CS_acceptor_porcent = prop.table(table(CS_acceptor$deltaG)) * 100
dados_hist_CS_acceptor = data.frame(delta_g = as.numeric(names(CS_acceptor_porcent)), CS_acceptor_porcent)

# ALTA (acceptor)
ALTA_porcent = prop.table(table(ALTA_acceptor$deltaG)) * 100
dados_hist_ALTA = data.frame(delta_g = as.numeric(names(ALTA_porcent)), ALTA_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Constitutive Acceptor"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Stability distributions - Constitutive Acceptor vs ALTA",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Constitutive Acceptor" = "royalblue", "ALTA" = "pink"))

# INT (acceptor)
INT_acceptor_porcent = prop.table(table(INT_acceptor$deltaG)) * 100
dados_hist_INT_acceptor = data.frame(delta_g = as.numeric(names(INT_acceptor_porcent)), INT_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Constitutive Acceptor"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_acceptor, aes(x = delta_g, y = INT_acceptor_porcent, fill = "INT acceptor"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Stability distributions - Constitutive Acceptor vs INT Acceptor",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Constitutive Acceptor" = "royalblue", "INT acceptor" = "pink"))

# EX events (acceptor)
EX_acceptor_porcent = prop.table(table(EX_acceptor$deltaG)) * 100
dados_hist_EX_acceptor = data.frame(delta_g = as.numeric(names(EX_acceptor_porcent)), EX_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Constitutive Acceptor"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX acceptor"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Stability distributions - Constitutive Acceptor vs EX Acceptor",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Constitutive Acceptor" = "royalblue", "EX acceptor" = "pink"))





