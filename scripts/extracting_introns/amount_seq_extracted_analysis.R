library(ggplot2)

# Dados de exemplo
genomas <- c("Arabidopsis", "Human hg18", "Human hg38", "Maize")
CS <- c(29886,9043, 128437, 4410)
ALTA <- c(85049,3634,11729,24748)
ALTD <- c(43643,5332,19849,17871)
INT <- c(115163,5150,34905,46416)
EX <- c(92764,24511,85268,14531)

# Criar um dataframe com os dados
dados <- data.frame(Genoma = rep(genomas, 5),
                    Event = rep(c( "CS", "ALTA", "ALTD", "INT", "EX"), each = 4),
                    Valor = c(CS, ALTA, ALTD, INT, EX))

cores = c("lightpink", "royalblue", "skyblue1", "khaki", "palegreen2")

# Criar o gráfico de barras usando ggplot2
grafico = ggplot(dados, aes(x = Genoma, y = Valor, fill = Event)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Genome", y = "Sequences amount") +
  scale_fill_manual(values = cores) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  theme_minimal() 

grafico + ggtitle("Sequences amount per genome to predict secondary structure") 
# + geom_text(aes(label = Valor), position = position_stack(vjust = 0.5), color = "white", size = 3, vjust = -0.2)

