
library(ggplot2)
CS_ALTD <- c(39.87, 41.83, 42.17)
dCS_INT <- c(48.82, 49.67, 49.15)
dCS_EX <- c(48.25, 48.97, 50.13)
dAlternative_Constitutive <- c(47.77,47.55,48.04)

CS_ALTA <- c(40.83, 40.72, 41.33)
CS_INT <- c(48.17, 49, 47.95)
CS_EX <- c(49.43, 50.27, 48.67)
aAlternative_Constitutive <- c(48.07, 47.02, 48.02)


# ... e assim por diante para os outros conjuntos
conjuntos <- data.frame(
  Conjunto = rep(c("CS vs ALTD", "CS vs INT (donor)", "CS vs EX (donor)", "CS x ALTERNATIVE (donor)", "CS vs ALTA", "CS vs INT (acceptor)", "CS vs EX (acceptor)", "CS x ALTERNATIVE (acceptor)"), each = 3),
  seed = rep(c("123", "231", "312"), times = 8),
  Dados = c(CS_ALTD,dCS_INT, dCS_EX, dAlternative_Constitutive, CS_ALTA, CS_INT, CS_EX, aAlternative_Constitutive)
)


p<-ggplot(conjuntos, aes(x=Conjunto, y=Dados, fill=seed)) +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.5), dotsize=0.9) +
  ylim(0,100) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Categories", y = "OOB error", title = "A. thaliana") 
p

p+scale_fill_manual(values=c("lightpink", "palegreen2", "skyblue2"))


CS_ALTD <- c(36.9, 38.57, 37.93)
dCS_INT <- c(36.7, 37.15, 34.98)
dCS_EX <- c(38.3, 39.18, 37.18)
dAlternative_Constitutive <- c(40.3, 40.18, 40.22)

CS_ALTA <- c(42.37, 42.37, 41.68)
CS_INT <- c(37.27, 37.92, 38)
CS_EX <- c(37.35, 38.83, 38.72)
aAlternative_Constitutive <- c(41.25, 42, 42.47)


# ... e assim por diante para os outros conjuntos
conjuntos <- data.frame(
  Conjunto = rep(c("CS vs ALTD", "CS vs INT (donor)", "CS vs EX (donor)", "CS x ALTERNATIVE (donor)", "CS vs ALTA", "CS vs INT (acceptor)", "CS vs EX (acceptor)", "CS x ALTERNATIVE (acceptor)"), each = 3),
  seed = rep(c("123", "231", "312"), times = 8),
  Dados = c(CS_ALTD,dCS_INT, dCS_EX, dAlternative_Constitutive, CS_ALTA, CS_INT, CS_EX, aAlternative_Constitutive)
)


p<-ggplot(conjuntos, aes(x=Conjunto, y=Dados, fill=seed)) +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.5), dotsize=0.9) +
  ylim(0,100) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Categories", y = "OOB error", title = "Maize") 
p

p+scale_fill_manual(values=c("lightpink", "palegreen2", "skyblue2"))




CS_ALTD <- c(30.25,30.32,29.85)
dCS_INT <- c(37.08,35.38,35.03)
dCS_EX <- c(42.32,41.53,42.4)
dAlternative_Constitutive <- c(39.3,38.78,37.92)

CS_ALTA <- c(27.97,26.82,26.7)
CS_INT <- c(35.22,36.57,34.65)
CS_EX <- c(42.17,40.88,40.95)
aAlternative_Constitutive <- c(37.88,38.75,38.68)


# ... e assim por diante para os outros conjuntos
conjuntos <- data.frame(
  Conjunto = rep(c("CS vs ALTD", "CS vs INT (donor)", "CS vs EX (donor)", "CS x ALTERNATIVE (donor)", "CS vs ALTA", "CS vs INT (acceptor)", "CS vs EX (acceptor)", "CS x ALTERNATIVE (acceptor)"), each = 3),
  seed = rep(c("123", "231", "312"), times = 8),
  Dados = c(CS_ALTD,dCS_INT, dCS_EX, dAlternative_Constitutive, CS_ALTA, CS_INT, CS_EX, aAlternative_Constitutive)
)


p<-ggplot(conjuntos, aes(x=Conjunto, y=Dados, fill=seed)) +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.5), dotsize=0.9) +
  ylim(0,100) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Categories", y = "OOB error", title = "Human hg38") 
p

p+scale_fill_manual(values=c("lightpink", "palegreen2", "skyblue2"))
