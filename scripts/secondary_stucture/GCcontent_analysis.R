library(seqinr)
library(LncFinder)
library(ggplot2)

f <- read.fasta('/home/bia/sugarcane_introns_local/data/maize/generated/sequences/ALL_2017v31_maize_CS_acceptor.fa')
gc <- compute_GC(f)
g <- read.fasta('/home/bia/sugarcane_introns_local/data/maize/generated/sequences/ALL_2017v31_maize_EXsk_acceptor.fa')
gc2 <- compute_GC(g)

G1 = "Constitutivo"
G2 = "EX"
ylabel = 'Conteúdo de GC'
titulo = ""

ggplot() +
  geom_boxplot(data = gc, aes(x = G1, y = GC.content, fill = "Group 1")) +
  geom_boxplot(data = gc2, aes(x = G2, y = GC.content, fill = "Group 2")) +
  ylab("Conteúdo de GC") +
  xlab("Evento de Splicing") +
  ggtitle(titulo) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Group 1" = "pink", "Group 2" = "salmon"))

##########

library(seqinr)
library(LncFinder)
library(ggplot2)

f <- read.fasta('/home/bia/sugarcane_introns_local/data/maize/generated/sequences/ALL_2017v31_maize_CS_donor.fa')
gc <- compute_GC(f)
g <- read.fasta('/home/bia/sugarcane_introns_local/data/maize/generated/sequences/ALL_2017v31_maize_INT_donor.fa')
gc2 <- compute_GC(g)
h <- read.fasta('/home/bia/sugarcane_introns_local/data/maize/generated/sequences/ALL_2017v31_maize_EXsk_donor.fa')
gc3 <- compute_GC(h)
i <- read.fasta('/home/bia/sugarcane_introns_local/data/maize/generated/sequences/ALL_2017v31_maize_ALTD_donor.fa')
gc4 <- compute_GC(i)

G1 = "Constitutivo"
G2 = "INT"
G3 = "EX"
G4 = "ALTD"
ylabel = 'Conteúdo de GC'
titulo = ""

ggplot() +
  geom_boxplot(data = gc, aes(x = G1, y = GC.content, fill = "Group 1")) +
  geom_boxplot(data = gc2, aes(x = G2, y = GC.content, fill = "Group 2")) +
  geom_boxplot(data = gc3, aes(x = G3, y = GC.content, fill = "Group 3")) +
  geom_boxplot(data = gc4, aes(x = G4, y = GC.content, fill = "Group 4")) +
  ylab("Conteúdo de GC") +
  xlab("Evento de Splicing") +
  ggtitle(titulo) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Group 1" = "lightblue", "Group 2" = "royalblue", "Group 3" = "darkblue", "Group 4" = "skyblue")) +
  ylim(0,1)




library(seqinr)
library(LncFinder)
library(ggplot2)

f <- read.fasta('/home/bia/sugarcane_introns_local/data/Human_genome/generated/sequences/ALL_hg18_human_CS_donor.fa')
gc <- compute_GC(f)
g <- read.fasta('/home/bia/sugarcane_introns_local/data/Human_genome//generated/sequences/ALL_hg18_human_INT_donor.fa')
gc2 <- compute_GC(g)
h <- read.fasta('/home/bia/sugarcane_introns_local/data/Human_genome/generated/sequences/ALL_hg18_human_EXsk_donor.fa')
gc3 <- compute_GC(h)

G1 = "Constitutivo"
G2 = "INT"
G3 = "EX"
ylabel = 'Conteúdo de GC'
titulo = ""

ggplot() +
  geom_boxplot(data = gc, aes(x = G1, y = GC.content, fill = "Group 1")) +
  geom_boxplot(data = gc2, aes(x = G2, y = GC.content, fill = "Group 2")) +
  geom_boxplot(data = gc3, aes(x = G3, y = GC.content, fill = "Group 3")) +
  ylab("Conteúdo de GC") +
  xlab("Evento de Splicing") +
  ggtitle(titulo) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("Group 1" = "lightblue", "Group 2" = "royalblue", "Group 3" = "darkblue")) +
  ylim(0,1)







