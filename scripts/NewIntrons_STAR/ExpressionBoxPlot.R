library(ggplot2)
CTBE = read.csv("/home/bia/sugarcane_introns_local/NewIntrons_STAR/CTBE_NewIntrons_STAR.csv", header=TRUE)
Souza = read.csv("/home/bia/sugarcane_introns_local/NewIntrons_STAR/Souza_NewIntrons_STAR.csv", header=TRUE)

CTBE$Genoma <- "CTBE"
Souza$Genoma <- "Souza"
combined_data <- rbind(CTBE, Souza)

ggplot(combined_data, aes(x = Genoma, y = mapped_sum, fill = Genoma)) +
  geom_boxplot() +
  labs(x = "Genoma", y = "Quantidade de Leituras Mapeadas") +
  geom_point(stat = "summary", fun = "mean", color = "tomato", size = 3) +
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = c("CTBE" = "skyblue2", "Souza" = "lightgreen")) +
  ggtitle("Quantidade de leituras nos novos introns")
