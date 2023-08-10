
aCTBE = read.csv('/home/bia/sugarcane_introns_local/data/sugarcane/CTBE_struc_acceptor.csv', header = TRUE)
dCTBE = read.csv('/home/bia/sugarcane_introns_local/data/sugarcane/CTBE_struc_donor.csv', header = TRUE)
aSouza = read.csv('/home/bia/sugarcane_introns_local/data/sugarcane/Souza_struc_acceptor.csv', header = TRUE)
dSouza = read.csv('/home/bia/sugarcane_introns_local/data/sugarcane/Souza_struc_donor.csv', header = TRUE)

aCTBE_porcent = prop.table(table(aCTBE$deltaG)) * 100
dados_hist_aCTBE = data.frame(delta_g = as.numeric(names(aCTBE_porcent)), aCTBE_porcent)

ggplot() +
  geom_bar(data = dados_hist_aCTBE, aes(x = delta_g, y = aCTBE_porcent, fill = "acceptor"), stat = "identity", alpha = 0.5,color = "royalblue") +
  labs(title = "Sugarcane deltaG distribution - Acceptor - CTBE",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("acceptor" = "royalblue"))

aSouza_porcent = prop.table(table(aSouza$deltaG)) * 100
dados_hist_aSouza = data.frame(delta_g = as.numeric(names(aSouza_porcent)), aSouza_porcent)
ggplot() +
  geom_bar(data = dados_hist_aSouza, aes(x = delta_g, y = aSouza_porcent, fill = "acceptor"), stat = "identity", alpha = 0.5,color = "lightgreen") +
  labs(title = "Sugarcane deltaG distribution - Acceptor - Souza",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("acceptor" = "lightgreen"))


dCTBE_porcent = prop.table(table(dCTBE$deltaG)) * 100
dados_hist_dCTBE = data.frame(delta_g = as.numeric(names(dCTBE_porcent)), dCTBE_porcent)

ggplot() +
  geom_bar(data = dados_hist_dCTBE, aes(x = delta_g, y = dCTBE_porcent, fill = "donor"), stat = "identity", alpha = 0.5,color = "royalblue") +
  labs(title = "Sugarcane deltaG distribution - donor - CTBE",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("donor" = "royalblue"))

dSouza_porcent = prop.table(table(dSouza$deltaG)) * 100
dados_hist_dSouza = data.frame(delta_g = as.numeric(names(dSouza_porcent)), dSouza_porcent)
ggplot() +
  geom_bar(data = dados_hist_dSouza, aes(x = delta_g, y = dSouza_porcent, fill = "donor"), stat = "identity", alpha = 0.5,color = "lightgreen") +
  labs(title = "Sugarcane deltaG distribution - Donor - Souza",
       x = "Valor do (kcal/mol)",
       y = "Porcentagem de Ocorrência",
       fill = "Evento") +
  xlim(-100, 5) +
  ylim(0,0.5) +
  scale_fill_manual(values = c("donor" = "lightgreen"))


## Coeficiente de co relação 
library(ggpubr)
CTBE = read.csv('/home/bia/sugarcane_introns_local/data/sugarcane/CTBE.final_table.csv', header = TRUE)
Souza = read.csv('/home/bia/sugarcane_introns_local/data/sugarcane/Souza.final_table.csv', header = TRUE)

CTBE$acceptor_deltaG = aCTBE$deltaG
CTBE$donor_deltaG = dCTBE$deltaG

Souza$acceptor_deltaG = aSouza$deltaG
Souza$donor_deltaG = dSouza$deltaG

ggscatter(CTBE, x=acceptor_deltaG, y=, add = "reg.line", conf.int=TRUE,
          cor.coef)

