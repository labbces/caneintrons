##############################################
#######---------- A THALIANA ----------#######
##############################################

# loading analysis csv files
dCS = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_CS_donor.csv', header = TRUE)
aCS = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_CS_acceptor.csv', header = TRUE)
ALTA = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_ALTA_acceptor.csv', header = TRUE)
ALTD = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_ALTD_donor.csv', header = TRUE)
dINT = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_INT_donor.csv', header = TRUE)
aINT = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_INT_acceptor.csv', header = TRUE)
dEX = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_EXsk_donor.csv', header = TRUE)
aEX = read.csv('/home/bia/sugarcane_introns_local/data/Athaliana_genome/generated/features_grapple_din/TAIR10_ALL_AtRTDv2_QUASI_EXsk_acceptor.csv', header = TRUE)

# checking values type - all should be numeric or integer, chr shouldn't be here, change it for factor bellow
str(dCS)
str(aCS)
str(ALTA)
str(ALTD)
str(dINT)
str(aINT)
str(dEX)
str(aEX)

###############################################################################
### IN CASE chr class IN DF, CHANGE IT TO FACTOR USING THE FOLLOWING LINES ###
#
dCS$label = as.factor(dCS$label) 
aCS$label = as.factor(aCS$label) 
ALTA$label = as.factor(ALTA$label) 
ALTD$label = as.factor(ALTD$label) 
dINT$label = as.factor(dINT$label) 
aINT$label = as.factor(aINT$label) 
dEX$label = as.factor(dEX$label) 
aEX$label = as.factor(aEX$label) 

###############################################################################

##############################################################################
### IN CASE non needed COLUMN IN DF, use the following code to exclude it ###
# 
dCS = subset(dCS, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
aCS = subset(aCS, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
ALTA = subset(ALTA, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
ALTD = subset(ALTD, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
dINT = subset(dINT, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
aINT = subset(aINT, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
dEX = subset(dEX, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
aEX = subset(aEX, select = -c(freq,freq1.2,freq2.2,freq1.3,freq2.3,freq3.3,substru,Stem_comp))
###################################################################################

#########################################################################################
### IN CASE OF nonNA NA VALUES IN DF, IT IS POSSIBLE USE THE FOLLOWING CODE TO FIX IT ###
#
# dCS = rfImpute(label ~ ., data = dCS, iter=6) # in theory, iter should be ok between 4 and 6
# aCS = rfImpute(label ~ ., data = aCS, iter=6)
# ALTA = rfImpute(label ~ ., data = ALTA, iter=6)
# ALTD = rfImpute(label ~ ., data = ALTD, iter=6)
# dINT = rfImpute(label ~ ., data = dINT, iter=6)
# aINT = rfImpute(label ~ ., data = aINT, iter=6)
# dEX = rfImpute(label ~ ., data = dEX, iter=6)
# aEX = rfImpute(label ~ ., data = aEX, iter=6)
#########################################################################################

####################################################################################
### IN CASE OF LONG DF, IT IS POSSIBLE TO MAKE A SAMPLE USING THE FOLLOWING CODE ###
sample_size = 3000 # Change for the sample size you want
set.seed(30)

library(dplyr)

dCS = sample_n(dCS, sample_size)
aCS = sample_n(aCS, sample_size)
ALTA = sample_n(ALTA, sample_size)
ALTD = sample_n(ALTD, sample_size)
dINT = sample_n(dINT, sample_size)
aINT = sample_n(aINT, sample_size)
dEX = sample_n(dEX, sample_size)
aEX = sample_n(aEX, sample_size)
####################################################################################

# binding analysis DF for random forest running
aCS_ALTA = rbind(aCS, ALTA)
dCS_ALTD = rbind(dCS, ALTD)
dCS_INT = rbind(dCS, dINT)
aCS_INT = rbind(aCS, aINT)
dCS_EX = rbind(dCS, dEX)
aCS_EX = rbind(aCS, aEX)




# (almost) all together - separated by ss type
CS_DONOR = rbind(dCS, ALTD, dINT, dEX)
CS_ACCEPTOR = rbind(aCS, ALTA, aINT, aEX)
CS_ACCEPTOR = rbind(aCS, aINT, aEX)

# como parece ser um evento muito mais relacionado com cada ss e nao com o conjunto
# faria sentido unir o donor e acceptor na an[alise?

# event 
CS = merge(aCS, dCS, by = 'id', all = TRUE)
INT = merge(aINT, dINT, by = 'id', all = TRUE)
EX = merge(aEX, dEX, by = 'id', all = TRUE)
CS_INT = rbind(CS, INT)
CS_EX = rbind(CS, EX)

# (again, almost all, ALTD and ALTA not here) all together
all = rbind(CS, INT, EX)

# Random Forest
library(randomForest)
# help(randomForest)

# From pair biding 
aCS_ALTA$label = as.factor(aCS_ALTA$label) 
aCS_ALTA = subset(aCS_ALTA, select = -c(id))
randomForest(label ~ .,data=aCS_ALTA, ntree=500) # ! wrong? why?
str(ALTA)


dCS_ALTD = subset(dCS_ALTD, select = -c(id))
randomForest(label ~ .,data=dCS_ALTD, ntree=500)

dCS_INT = subset(dCS_INT, select = -c(id))
randomForest(label ~ .,data=dCS_INT, ntree=500)  

aCS_INT = subset(aCS_INT, select = -c(id))
randomForest(label ~ .,data=aCS_INT, ntree=500)

dCS_EX = subset(dCS_EX, select = -c(id))
randomForest(label ~ .,data=dCS_EX, ntree=500)

aCS_EX = subset(aCS_EX, select = -c(id))
randomForest(label ~ .,data=aCS_EX, ntree=500) ##

# From separated by SS type 
CS_DONOR$label = as.factor(CS_DONOR$label) 
CS_DONOR = subset(CS_DONOR, select = -c(id))
randomForest(label ~ .,data=CS_DONOR, ntree=500)

str(CS_ACCEPTOR)
CS_ACCEPTOR$label = as.factor(CS_ACCEPTOR$label) 
CS_ACCEPTOR = subset(CS_ACCEPTOR, select = -c(id))
randomForest(label ~ .,data=CS_ACCEPTOR, ntree=500)





sample_size_CS = 3000 # Change for the sample size you want
sample_size = 1000
set.seed(30)

library(dplyr)

dCS = sample_n(dCS, sample_size_CS)
aCS = sample_n(aCS, sample_size_CS)
ALTA = sample_n(ALTA, sample_size)
ALTD = sample_n(ALTD, sample_size)
dINT = sample_n(dINT, sample_size)
aINT = sample_n(aINT, sample_size)
dEX = sample_n(dEX, sample_size)
aEX = sample_n(aEX, sample_size)

CS_DONOR = rbind(dCS, ALTD, dINT, dEX)
CS_ACCEPTOR = rbind(aCS, ALTA, aINT, aEX)


CS_DONOR = subset(CS_DONOR, select = -c(id))
library(dplyr)
CS_DONOR$label = as.character(CS_DONOR$label) 
CS_DONOR["label"][CS_DONOR["label"] == "ALTA"] <- "ALTERNATIVO"
CS_DONOR["label"][CS_DONOR["label"] == "ALTD"] <- "ALTERNATIVO"
CS_DONOR["label"][CS_DONOR["label"] == "INT"] <- "ALTERNATIVO"
CS_DONOR["label"][CS_DONOR["label"] == "EXsk"] <- "ALTERNATIVO"
CS_DONOR$label = as.factor(CS_DONOR$label) 
r = randomForest(label ~ .,data=CS_DONOR, ntree=500)
r
varImpPlot(r)


CS_ACCEPTOR = subset(CS_ACCEPTOR, select = -c(id))
CS_ACCEPTOR$label = as.character(CS_ACCEPTOR$label) 
CS_ACCEPTOR["label"][CS_ACCEPTOR["label"] == "ALTA"] <- "ALTERNATIVO"
CS_ACCEPTOR["label"][CS_ACCEPTOR["label"] == "ALTD"] <- "ALTERNATIVO"
CS_ACCEPTOR["label"][CS_ACCEPTOR["label"] == "INT"] <- "ALTERNATIVO"
CS_ACCEPTOR["label"][CS_ACCEPTOR["label"] == "EXsk"] <- "ALTERNATIVO"
CS_ACCEPTOR$label = as.factor(CS_ACCEPTOR$label) 
r = randomForest(label ~ .,data=CS_ACCEPTOR, ntree=500)
r
varImpPlot(r)



