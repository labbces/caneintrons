title = c("ALTA", "ALTD", "donor INT", "acceptor INT", "donor EX", "acceptor EX", "All donor", 'All acceptor')
human_oob = c(41.12, 45.67, 38.15, 39.83, 50.38, 49.83, 66.56, 66.17)
maize_oob = c(46.12, 44.72, 0.88, 43.12, 44.33, 44.27, 45.51, 68.42)
athali_oob = c(NaN, 47.18, 49.68, 49.3, 51.1, 0.62, 74.61, 34.39) # hey, last one is missing ALTA!

values = c(human_oob, maize_oob, athali_oob)
source = c("human", "human", "human", "human", "human", "human", "human", "human", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "Athaliana", "Athaliana", "Athaliana", "Athaliana", "Athaliana", "Athaliana", "Athaliana", "Athaliana")
type = c(title, title, title)

df = data.frame(values, source, type, stringsAsFactors = TRUE)
df

library(ggplot2)
p = ggplot(df, aes(x=type, y=values, fill=source, label = values)) + geom_dotplot(binaxis = 'y',
                                                              binpositions="bygroup", 
                                                              dotsize = 1,
                                                              stackdir='center') + ylim(0, 100)

p+scale_fill_manual(values=c("lightgreen", "chocolate1", "gold")) +
  ggtitle("Random Forest OOB error") + theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "OOB error", x = "alternative splicing event compared to constitutive splicing") +
  geom_text(size = 3, hjust = -0.5)
  