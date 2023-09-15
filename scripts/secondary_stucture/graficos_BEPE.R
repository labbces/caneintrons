
library(ggplot2)


################## HUMAN HG18 ################
CS_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_CS_donor.csv')
CS_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_CS_acceptor.csv')
ALTD_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_ALTD_donor.csv')
ALTA_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_ALTA_acceptor.csv')
INT_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_INT_donor.csv')
INT_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_INT_acceptor.csv')
EX_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_EXsk_donor.csv')
EX_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg18/features/ALL_hg18_human_EXsk_acceptor.csv')

#DonorHG18###############
CS_donor_porcent = prop.table(table(CS_donor$deltaG)) * 100
dados_hist_CS_donor = data.frame(delta_g = as.numeric(names(CS_donor_porcent)), CS_donor_porcent)

ALTD_porcent = prop.table(table(ALTD_donor$deltaG)) * 100
dados_hist_ALTD = data.frame(delta_g = as.numeric(names(ALTD_porcent)), ALTD_porcent)

ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador Constitutivo"), stat = "identity", alpha = 0.5,color = "royalblue") +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.5, color = "pink") +
  labs(title = "Humano HG18 - Doador Constitutivo vs ALTD",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Doador Constitutivo" = "royalblue", "ALTD" = "pink"))


EX_donor_porcent = prop.table(table(EX_donor$deltaG)) * 100
dados_hist_EX_donor = data.frame(delta_g = as.numeric(names(EX_donor_porcent)), EX_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador Constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX Doador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG18 - Doanor Constitutivo vs EX Doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Doador Constitutivo" = "royalblue", "EX Doador" = "pink"))

INT_donor_porcent = prop.table(table(INT_donor$deltaG)) * 100
dados_hist_INT_donor = data.frame(delta_g = as.numeric(names(INT_donor_porcent)), INT_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador Constitutivo"), stat = "identity", color = "royalblue") +
  geom_bar(data = dados_hist_INT_donor, aes(x = delta_g, y = INT_donor_porcent, fill = "INT Doador"), stat = "identity", color = "lightpink2") +
  labs(title = "Humano HG18 - Doador Constitutivo vs INT Doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Doador Constitutivo" = "royalblue", "INT Doador" = "lightpink2"))



ggplot() +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX Doador"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Humano HG18 - ALTD vs EX Doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("ALTD" = "darkseagreen3", "EX Doador" = "palevioletred1"))






CS_donor$origem = "CS"
ALTD_donor$origem = "ALTD"
EX_donor$origem = "EX"
INT_donor$origem = "INT"
combined_data <- rbind(CS_donor, ALTD_donor, INT_donor, EX_donor)
combined_data18 =  combined_data
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTD", "EX", "INT")) +
  ylim(-100, 0) +
  geom_hline(yintercept = median(CS_donor$deltaG), linetype = "dashed", color = "red") + 
  ggtitle("Human HG18 - Doador") + xlab("Evento") + ylab("kcal/mol") + stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")

median(CS_donor$deltaG)
median(ALTD_donor$deltaG)
median(EX_donor$deltaG)
median(INT_donor$deltaG)


#AcceptorHG18####
CS_acceptor_porcent = prop.table(table(CS_acceptor$deltaG)) * 100
dados_hist_CS_acceptor = data.frame(delta_g = as.numeric(names(CS_acceptor_porcent)), CS_acceptor_porcent)

ALTA_porcent = prop.table(table(ALTA_acceptor$deltaG)) * 100
dados_hist_ALTA = data.frame(delta_g = as.numeric(names(ALTA_porcent)), ALTA_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG18 - Aceitador constitutivo vs ALTA",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "ALTA" = "pink"))

INT_acceptor_porcent = prop.table(table(INT_acceptor$deltaG)) * 100
dados_hist_INT_acceptor = data.frame(delta_g = as.numeric(names(INT_acceptor_porcent)), INT_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_acceptor, aes(x = delta_g, y = INT_acceptor_porcent, fill = "INT Aceitador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG18 - Aceitador constitutivo vs INT Aceitador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "INT Aceitador" = "pink"))

EX_acceptor_porcent = prop.table(table(EX_acceptor$deltaG)) * 100
dados_hist_EX_acceptor = data.frame(delta_g = as.numeric(names(EX_acceptor_porcent)), EX_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX Aceitador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG18 - Aceitador constitutivo vs EX Aceitador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "EX Aceitador" = "pink"))

CS_acceptor$origem = "CS"
ALTA_acceptor$origem = "ALTA"
EX_acceptor$origem = "EX"
INT_acceptor$origem = "INT"
combined_data <- rbind(CS_acceptor, ALTA_acceptor, INT_acceptor, EX_acceptor)
combined_data18a = combined_data
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTA", "EX", "INT")) +
  ylim(-100, 0) + coord_flip() +
  geom_hline(yintercept = median(CS_acceptor$deltaG), linetype = "dashed", color = "red") + ggtitle("Human HG18 - Aceitador") + xlab("Evento") + ylab("kcal/mol") + stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")

ggplot() +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX Acceptor"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Humano HG18 - ALTA vs EX Acceptor",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("ALTA" = "darkseagreen3", "EX Acceptor" = "palevioletred1"))


################## HUMAN HG38 ################
CS_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_CS_donor.csv')
CS_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_CS_acceptor.csv')
ALTD_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_ALTD_donor.csv')
ALTA_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_ALTA_acceptor.csv')
INT_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_INT_donor.csv')
INT_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_INT_acceptor.csv')
EX_donor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_EXsk_donor.csv')
EX_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Human_genome/generated_hg38/ALL_hg38_human_EXsk_acceptor.csv')

#DonorHG38###################
CS_donor_porcent = prop.table(table(CS_donor$deltaG)) * 100
dados_hist_CS_donor = data.frame(delta_g = as.numeric(names(CS_donor_porcent)), CS_donor_porcent)

ALTD_porcent = prop.table(table(ALTD_donor$deltaG)) * 100
dados_hist_ALTD = data.frame(delta_g = as.numeric(names(ALTD_porcent)), ALTD_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG38 - Doador constitutivo vs ALTD",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "ALTD" = "pink"))

INT_donor_porcent = prop.table(table(INT_donor$deltaG)) * 100
dados_hist_INT_donor = data.frame(delta_g = as.numeric(names(INT_donor_porcent)), INT_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_donor, aes(x = delta_g, y = INT_donor_porcent, fill = "INT doador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG38 - Doador constitutivo vs INT doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "INT doador" = "pink"))

EX_donor_porcent = prop.table(table(EX_donor$deltaG)) * 100
dados_hist_EX_donor = data.frame(delta_g = as.numeric(names(EX_donor_porcent)), EX_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX doador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG38 - Doador constitutivo vs EX doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "EX doador" = "pink"))

CS_donor$origem = "CS"
ALTD_donor$origem = "ALTD"
EX_donor$origem = "EX"
INT_donor$origem = "INT"
combined_data <- rbind(CS_donor, ALTD_donor, INT_donor, EX_donor)
combined_data38 =  combined_data
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTD", "EX", "INT")) +
  ylim(-100, 0) + coord_flip() +
  geom_hline(yintercept = median(CS_donor$deltaG), linetype = "dashed", color = "red") + ggtitle("Human HG38 - Doador") + xlab("Evento") + ylab("kcal/mol") + stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")

median(CS_donor$deltaG)
median(ALTD_donor$deltaG)
median(EX_donor$deltaG)
median(INT_donor$deltaG)
median(CS_donor$deltaG)
median(ALTD_donor$deltaG)
median(EX_donor$deltaG)
median(INT_donor$deltaG)

#AcceptorHG38###########
CS_acceptor_porcent = prop.table(table(CS_acceptor$deltaG)) * 100
dados_hist_CS_acceptor = data.frame(delta_g = as.numeric(names(CS_acceptor_porcent)), CS_acceptor_porcent)


ALTA_porcent = prop.table(table(ALTA_acceptor$deltaG)) * 100
dados_hist_ALTA = data.frame(delta_g = as.numeric(names(ALTA_porcent)), ALTA_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG38 - Aceitador constitutivo vs ALTA",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "ALTA" = "pink"))


INT_acceptor_porcent = prop.table(table(INT_acceptor$deltaG)) * 100
dados_hist_INT_acceptor = data.frame(delta_g = as.numeric(names(INT_acceptor_porcent)), INT_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_acceptor, aes(x = delta_g, y = INT_acceptor_porcent, fill = "INT aceitador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG38 - Aceitador constitutivo vs INT aceitador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "INT aceitador" = "pink"))

EX_acceptor_porcent = prop.table(table(EX_acceptor$deltaG)) * 100
dados_hist_EX_acceptor = data.frame(delta_g = as.numeric(names(EX_acceptor_porcent)), EX_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX aceitador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Humano HG38 - Aceitador constitutivo vs EX Aceitador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "EX aceitador" = "pink"))


CS_acceptor$origem = "CS"
ALTA_acceptor$origem = "ALTA"
EX_acceptor$origem = "EX"
INT_acceptor$origem = "INT"
combined_data <- rbind(CS_acceptor, ALTA_acceptor, INT_acceptor, EX_acceptor)
combined_data38a = combined_data
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTA", "EX", "INT")) +
  ylim(-100, 0) + coord_flip() +
  geom_hline(yintercept = median(CS_acceptor$deltaG), linetype = "dashed", color = "red") + ggtitle("Human HG38 - Aceitador") + xlab("Evento") + ylab("kcal/mol") + stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")


ggplot() +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX Doador"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Humano HG38 - ALTD vs EX Doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("ALTD" = "darkseagreen3", "EX Doador" = "palevioletred1"))


ggplot() +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX Acceptor"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Humano HG38 - ALTA vs EX Acceptor",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("ALTA" = "darkseagreen3", "EX Acceptor" = "palevioletred1"))




#
#
#
library(dplyr)
hg18donor = combined_data18[, c("deltaG", "label")]
hg18donor$origem = "hg18"
hg38donor =  combined_data38[, c("deltaG", "label")]
hg38donor$origem = "hg38"

ave = rbind(hg18donor, hg38donor)
ggplot(ave, aes(x=label, y=deltaG, fill=origem)) + 
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTD", "EXsk", "INT")) +
  ylim(-100, 0) + 
  ggtitle("Doador") + xlab("Evento") + ylab("kcal/mol") + 
  scale_fill_manual(values= c("hg38"="#88d4e5", "hg18" = "#fdb827ff"))

hg18acceptor = combined_data18a[, c("deltaG", "label")]
hg18acceptor$origem = "hg18"
hg38acceptor =  combined_data38a[, c("deltaG", "label")]
hg38acceptor$origem = "hg38"

aves = rbind(hg18acceptor, hg38acceptor)
ggplot(aves, aes(x=label, y=deltaG, fill=origem)) + 
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTA", "EXsk", "INT")) +
  ylim(-100, 0) + 
  ggtitle("Aceitador") + xlab("Evento") + ylab("kcal/mol") + 
  scale_fill_manual(values= c("hg38"="#88d4e5", "hg18" = "#fdb827ff"))

  #
#
################## ARABIDOPSIS ################

CS_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_CS_donor.csv')
CS_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_CS_acceptor.csv')
ALTD_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_ALTD_donor.csv')
ALTA_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_ALTA_acceptor.csv')
INT_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_INT_donor.csv')
INT_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_INT_acceptor.csv')
EX_donor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_EXsk_donor.csv')
EX_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features/TAIR10_ALL_AtRTDv2_QUASI_EXsk_acceptor.csv')

#Donor#################
CS_donor_porcent = prop.table(table(CS_donor$deltaG)) * 100
dados_hist_CS_donor = data.frame(delta_g = as.numeric(names(CS_donor_porcent)), CS_donor_porcent)

ALTD_porcent = prop.table(table(ALTD_donor$deltaG)) * 100
dados_hist_ALTD = data.frame(delta_g = as.numeric(names(ALTD_porcent)), ALTD_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "A. thaliana - Doador constitutivo vs ALTD",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "ALTD" = "pink"))

# INT doador
INT_donor_porcent = prop.table(table(INT_donor$deltaG)) * 100
dados_hist_INT_donor = data.frame(delta_g = as.numeric(names(INT_donor_porcent)), INT_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_donor, aes(x = delta_g, y = INT_donor_porcent, fill = "INT doador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "A. thaliana - Doador constitutivo vs INT doador",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "INT doador" = "pink"))

# EX doador
EX_donor_porcent = prop.table(table(EX_donor$deltaG)) * 100
dados_hist_EX_donor = data.frame(delta_g = as.numeric(names(EX_donor_porcent)), EX_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX doador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "A. thaliana - Doador constitutivo vs EX doador",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "EX doador" = "pink"))


CS_donor$origem = "CS"
ALTD_donor$origem = "ALTD"
EX_donor$origem = "EX"
INT_donor$origem = "INT"
combined_data <- rbind(CS_donor, ALTD_donor, INT_donor, EX_donor)
maized = combined_data
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTD", "EX", "INT")) +
  ylim(-100, 0) + coord_flip() +
   ggtitle("Arabidopsis - Doador") + xlab("Evento") + ylab("kcal/mol") +
  scale_fill_manual(values= c("CS"="#88d4e5", "EX" = "#fdb827ff", "ALTD" = "#fdb827ff",  "INT" = "#fdb827ff"))


#Aceptor###########
CS_acceptor_porcent = prop.table(table(CS_acceptor$deltaG)) * 100
dados_hist_CS_acceptor = data.frame(delta_g = as.numeric(names(CS_acceptor_porcent)), CS_acceptor_porcent)

# ALTA (acceptor)
ALTA_porcent = prop.table(table(ALTA_acceptor$deltaG)) * 100
dados_hist_ALTA = data.frame(delta_g = as.numeric(names(ALTA_porcent)), ALTA_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "A. thaliana - Aceitador constitutivo vs ALTA",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "ALTA" = "pink"))

# INT (acceptor)
INT_acceptor_porcent = prop.table(table(INT_acceptor$deltaG)) * 100
dados_hist_INT_acceptor = data.frame(delta_g = as.numeric(names(INT_acceptor_porcent)), INT_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_acceptor, aes(x = delta_g, y = INT_acceptor_porcent, fill = "INT"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "A. thaliana - Aceitador constitutivo vs INT",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "INT" = "pink"))

# EX events (acceptor)
EX_acceptor_porcent = prop.table(table(EX_acceptor$deltaG)) * 100
dados_hist_EX_acceptor = data.frame(delta_g = as.numeric(names(EX_acceptor_porcent)), EX_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX aceitador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "A. thaliana - Aceitador constitutivo vs EX aceitador",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "EX aceitador" = "pink"))

CS_acceptor$origem = "CS"
ALTA_acceptor$origem = "ALTA"
EX_acceptor$origem = "EX"
INT_acceptor$origem = "INT"
combined_data <- rbind(CS_acceptor, ALTA_acceptor, INT_acceptor, EX_acceptor)
maizea = combined_data
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTA", "EX", "INT")) +
  ylim(-100, 0) + coord_flip() +
  ggtitle("Arabidopsis - Aceitador") + xlab("Evento") + ylab("kcal/mol") +
  scale_fill_manual(values= c("CS"="#88d4e5", "EX" = "#fdb827ff", "ALTA" = "#fdb827ff",  "INT" = "#fdb827ff"))


ggplot() +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX Doador"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Arabidopsis - ALTD vs EX Doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("ALTD" = "darkseagreen3", "EX Doador" = "palevioletred1"))


ggplot() +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX Acceptor"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Arabidopsis - ALTA vs EX Acceptor",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("ALTA" = "darkseagreen3", "EX Acceptor" = "palevioletred1"))




#
#
#
################## Milho ################
CS_donor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_CS_donor.csv')
CS_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_CS_acceptor.csv')
ALTD_donor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_ALTD_donor.csv')
ALTA_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_ALTA_acceptor.csv')
INT_donor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_INT_donor.csv')
INT_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_INT_acceptor.csv')
EX_donor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_EXsk_donor.csv')
EX_acceptor = read.csv('/home/bia/sugarcane_introns_local/data/maize/generated/features_grapple/ALL_2017v31_maize_EXsk_acceptor.csv')

#Donor################
CS_donor_porcent = prop.table(table(CS_donor$deltaG)) * 100
dados_hist_CS_donor = data.frame(delta_g = as.numeric(names(CS_donor_porcent)), CS_donor_porcent)

# ALTD
ALTD_porcent = prop.table(table(ALTD_donor$deltaG)) * 100
dados_hist_ALTD = data.frame(delta_g = as.numeric(names(ALTD_porcent)), ALTD_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Milho - Doador constitutivo vs ALTD",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "ALTD" = "pink"))

# INT donor
INT_donor_porcent = prop.table(table(INT_donor$deltaG)) * 100
dados_hist_INT_donor = data.frame(delta_g = as.numeric(names(INT_donor_porcent)), INT_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_donor, aes(x = delta_g, y = INT_donor_porcent, fill = "INT donor"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Milho - Doador constitutivo vs INT donor",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "INT donor" = "pink"))

# EX donor
median(dados_hist_EX_donor$delta_g)
median(dados_hist_INT_donor$delta_g)
median(dados_hist_CS_donor$delta_g)
EX_donor_porcent = prop.table(table(EX_donor$deltaG)) * 100
dados_hist_EX_donor = data.frame(delta_g = as.numeric(names(EX_donor_porcent)), EX_donor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_donor, aes(x = delta_g, y = CS_donor_porcent, fill = "Doador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX donor"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Milho - Doador constitutivo vs EX donor",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Doador constitutivo" = "royalblue", "EX donor" = "pink"))

CS_donor$origem = "CS"
ALTD_donor$origem = "ALTD"
EX_donor$origem = "EX"
INT_donor$origem = "INT"
combined_data <- rbind(CS_donor, ALTD_donor, INT_donor, EX_donor)
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTD", "EX", "INT")) +
  ylim(-100, 0) + coord_flip() +
  geom_hline(yintercept = median(CS_donor$deltaG), linetype = "dashed", color = "red") + ggtitle("Milho -  Doador") + xlab("Evento") + ylab("kcal/mol") + stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")

#Acceptor####################
CS_acceptor_porcent = prop.table(table(CS_acceptor$deltaG)) * 100
dados_hist_CS_acceptor = data.frame(delta_g = as.numeric(names(CS_acceptor_porcent)), CS_acceptor_porcent)

# ALTA (acceptor)
ALTA_porcent = prop.table(table(ALTA_acceptor$deltaG)) * 100
dados_hist_ALTA = data.frame(delta_g = as.numeric(names(ALTA_porcent)), ALTA_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Milho - Aceitador constitutivo vs ALTA",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "ALTA" = "pink"))

# INT (acceptor)
INT_acceptor_porcent = prop.table(table(INT_acceptor$deltaG)) * 100
dados_hist_INT_acceptor = data.frame(delta_g = as.numeric(names(INT_acceptor_porcent)), INT_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_INT_acceptor, aes(x = delta_g, y = INT_acceptor_porcent, fill = "INT aceitador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Milho - Aceitador constitutivo vs INT aceitador",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "INT aceitador" = "pink"))

# EX events (acceptor)
EX_acceptor_porcent = prop.table(table(EX_acceptor$deltaG)) * 100
dados_hist_EX_acceptor = data.frame(delta_g = as.numeric(names(EX_acceptor_porcent)), EX_acceptor_porcent)
ggplot() +
  geom_bar(data = dados_hist_CS_acceptor, aes(x = delta_g, y = CS_acceptor_porcent, fill = "Aceitador constitutivo"), stat = "identity", alpha = 1, color = "royalblue") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX aceitador"), stat = "identity", alpha = 0.1, color = "pink") +
  labs(title = "Milho - Aceitador constitutivo vs EX aceitador",
       x = "Valor do Delta G",
       y = "Porcentagem de Ocorrência",
       fill = "Event") +
  xlim(-100, 5) +
  ylim(0,1) +
  scale_fill_manual(values = c("Aceitador constitutivo" = "royalblue", "EX aceitador" = "pink"))

CS_acceptor$origem = "CS"
ALTA_acceptor$origem = "ALTA"
EX_acceptor$origem = "EX"
INT_acceptor$origem = "INT"
combined_data <- rbind(CS_acceptor, ALTA_acceptor, INT_acceptor, EX_acceptor)
ggplot(data = combined_data, aes(x = origem, y = deltaG, fill=origem)) +
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTA", "EX", "INT")) +
  ylim(-100, 0) + coord_flip() +
  geom_hline(yintercept = median(CS_acceptor$deltaG), linetype = "dashed", color = "red") + ggtitle("Milho -  Aceitador") + xlab("Evento") + ylab("kcal/mol") + stat_summary(fun.y=mean, geom="point", shape=20, size=2, color="red", fill="red")
#
ggplot() +
  geom_bar(data = dados_hist_ALTD, aes(x = delta_g, y = ALTD_porcent, fill = "ALTD"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_donor, aes(x = delta_g, y = EX_donor_porcent, fill = "EX Doador"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Maize - ALTD vs EX Doador",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.8) +
  scale_fill_manual(values = c("ALTD" = "darkseagreen3", "EX Doador" = "palevioletred1"))


ggplot() +
  geom_bar(data = dados_hist_ALTA, aes(x = delta_g, y = ALTA_porcent, fill = "ALTA"), stat = "identity", alpha = 0.5, color = "darkseagreen3") +
  geom_bar(data = dados_hist_EX_acceptor, aes(x = delta_g, y = EX_acceptor_porcent, fill = "EX Acceptor"), stat = "identity", alpha = 0.1, color = "palevioletred1") +
  labs(title = "Maize - ALTA vs EX Acceptor",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.8) +
  scale_fill_manual(values = c("ALTA" = "darkseagreen3", "EX Acceptor" = "palevioletred1"))


