### Different temperatures 

########ATHALIANA######
##############LoadingCSV
# 37 
dCS = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_CS_donor.csv')
aCS = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_CS_acceptor.csv')
ALTD = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_ALTD_donor.csv')
ALTA = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_ALTA_acceptor.csv')
dINT = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_INT_donor.csv')
aINT = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_INT_acceptor.csv')
dEX = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_EXsk_donor.csv')
aEX = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_EXsk_acceptor.csv')

# 25
dCS25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_CS_donor.csv')
aCS25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_CS_acceptor.csv')
ALTD25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_ALTD_donor.csv')
ALTA25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_ALTA_acceptor.csv')
dINT25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_INT_donor.csv')
aINT25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_INT_acceptor.csv')
dEX25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_EXsk_donor.csv')
aEX25 = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_temp/TAIR10_ALL_AtRTDv2_QUASI_EXsk_acceptor.csv')

# # Comparing deltaG
library(ggplot2)
library(dplyr)

dCS$origem = "CS"
ALTD$origem = "ALTD"
dEX$origem = "EX"
dINT$origem = "INT"
combined_data <- rbind(dCS, ALTD, dEX, dINT)
donor37 = combined_data[, c("deltaG", "label")]
donor37$origem = "37º"

dCS25$origem = "CS"
ALTD25$origem = "ALTD"
dEX25$origem = "EX"
dINT25$origem = "INT"
combined_data25 <- rbind(dCS, ALTD, dEX, dINT)
donor25 = combined_data25[, c("deltaG", "label")]
donor25$origem = "25º"

ave = rbind(donor25, donor37)
ggplot(ave, aes(x=label, y=deltaG, fill=origem)) + 
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTD", "EXsk", "INT")) +
  ylim(-100, 0) + 
  ggtitle("Athaliana - Temperature comparison - Donor") + xlab("Evento") + ylab("kcal/mol") + 
  scale_fill_manual(values= c("25º"="#88d4e5", "37º" = "#fdb827ff"))



aCS$origem = "CS"
ALTA$origem = "ALTA"
aEX$origem = "EX"
aINT$origem = "INT"
combined_data <- rbind(aCS, ALTA, aEX, aINT)
acceptor37 = combined_data[, c("deltaG", "label")]
acceptor37$origem = "37º"

aCS25$origem = "CS"
ALTA25$origem = "ALTA"
aEX25$origem = "EX"
aINT25$origem = "INT"
combined_data25 <- rbind(aCS, ALTA, aEX, aINT)
acceptor25 = combined_data25[, c("deltaG", "label")]
acceptor25$origem = "25º"

ave = rbind(acceptor25, acceptor37)
ggplot(ave, aes(x=label, y=deltaG, fill=origem)) + 
  geom_boxplot() +
  scale_x_discrete(limits = c("CS", "ALTA", "EXsk", "INT")) +
  ylim(-100, 0) + 
  ggtitle("Athaliana - Temperature comparison - Acceptor") + xlab("Evento") + ylab("kcal/mol") + 
  scale_fill_manual(values= c("25º"="#88d4e5", "37º" = "#fdb827ff"))

