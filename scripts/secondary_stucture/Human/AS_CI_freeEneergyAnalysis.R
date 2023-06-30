library(ggplot2)

# How data files were generated:
# 1. Extract SS flanking sequences with extractor 
# 2. Run RNAfold secondary structure prediction 
# 3. Run RNAfold analysis to get seconary stucture features - deltaG.
# - Run for retained data only, constitutive data only and both labeled data


### ANALYSIS OF HUMAN DATA ###

# Loading data - retained and constitutive secondary structure databases
# CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_constitutive_Acceptor.csv')
# CI_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_constitutive_Donor.csv')
CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_cons_ALL_constitutive_Acceptor_intron.csv')
CI_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_cons_ALL_constitutive_Donor_intron.csv')

altFive_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altFivePrime_Acceptor_intron.csv')
altFive_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altFivePrime_Donor_intron.csv')
altThree_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altThreePrime_Acceptor_intron.csv')
altThree_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altThreePrime_Donor_intron.csv')

retainedIntron_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_retainedIntron_Acceptor_intron.csv')
retainedIntron_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_retainedIntron_Donor_intron.csv')
cassette_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_cassetteExon_Acceptor_intron.csv')
cassette_donor=read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_cassetteExon_Donor_intron.csv')

### Generating histograms ###
# altFive - Donor
ggplot(altFive_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 5' Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# AS - Acceptor
ggplot(altFive_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 5' Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# altThree - Donor
ggplot(altThree_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 3' Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# altThree - Acceptor
ggplot(altThree_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 3' Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# Intron retention
ggplot(retainedIntron_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Retained intron Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
ggplot(retainedIntron_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Retained intron Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# Constitutive - Donor
ggplot(CI_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") +
  xlim(-100, 5)
  # Constitutive - Acceptor
ggplot(CI_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

### Balanced data AS and CI ###
CI_balenced_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_constitutive_Donor.csv')
CI_balenced_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_constitutive_Acceptor.csv')

CI_balenced_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_cons_ALL_balanced_constitutive_Donor_intron.csv')
CI_balenced_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_cons_ALL_balanced_constitutive_Acceptor_intron.csv')

# alternative Five
CI_altFive_constitutive = rbind(CI_balenced_acceptor, altFive_donor)
ggplot(CI_altFive_constitutive, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Acceptor vs Alt5 correct)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_alternative = rbind(CI_balenced_acceptor, altFive_acceptor)
ggplot(CI_altFive_alternative, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Acceptor vs Alt5 alternative)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_acceptor = rbind(CI_balenced_acceptor, altFive_acceptor)
ggplot(CI_altFive_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative 5 x Constitutive donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_acceptor = rbind(CI_balenced_acceptor, altFive_donor)
ggplot(CI_altFive_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative 5 - Donor splice site",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_donor = rbind(CI_balenced_acceptor, altFive_acceptor)
ggplot(CI_altFive_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# alternative Three
CI_altThree_constitutive = rbind(CI_balenced_donor, altThree_acceptor)
ggplot(CI_altThree_constitutive, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Donor vs Alt3 correct)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_alternative = rbind(CI_balenced_donor, altThree_donor)
ggplot(CI_altThree_alternative, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative (CS Donor vs Alt3 alternative)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_donor = rbind(CI_balenced_acceptor, altThree_donor)
ggplot(CI_altThree_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative 3' - Acceptor Splice site",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_acceptor = rbind(CI_balenced_donor, altThree_acceptor)
ggplot(CI_altThree_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - CS donor x AltThree acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_acceptor = rbind(CI_balenced_donor, altThree_donor)
ggplot(CI_altThree_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - CS donor x AltThree donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# IR
CI_IR_acceptor = rbind(retainedIntron_acceptor, CI_balenced_acceptor)
ggplot(CI_IR_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_IR_donor = rbind(CI_balenced_donor, retainedIntron_donor)
ggplot(CI_IR_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# Cassette
cassette_acceptor_balanced = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_cassetteExon_Acceptor_intron.csv')
cassette_donor_balanced = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_cassetteExon_Donor_intron.csv')
CI_cassette_aceptor = rbind(CI_balenced_acceptor, cassette_donor_balanced)
ggplot(CI_cassette_aceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exons - Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_cassette_donor = rbind(CI_balenced_donor, cassette_acceptor_balanced)
ggplot(CI_cassette_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exon - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

### All together ### menos o cssette
# Acceptor
all_together_acceptor= rbind(CI_balenced_acceptor, altThree_acceptor, altFive_acceptor, retainedIntron_acceptor )
ggplot(all_together_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.3) +
  labs(title="Stability distributions - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# Donor
all_together_donor= rbind(CI_balenced_donor, altThree_donor, altFive_donor, retainedIntron_donor)
ggplot(all_together_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.3) +
  labs(title="Stability distributions - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")


### ANALYSIS OF A. THALIANA DATA ###

# Loading data - retained and constitutive secondary structure databases
CI_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/RNAfold/TAIR10_cons_ALL_constitutive_Acceptor_intron.RNAfold')
CI_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Athaliana_genome/RNAfold/TAIR10_cons_ALL_constitutive_Donor_intron.RNAfold')

altFive_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altFivePrime_Acceptor_intron.csv')
altFive_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altFivePrime_Donor_intron.csv')
altThree_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altThreePrime_Acceptor_intron.csv')
altThree_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_altThreePrime_Donor_intron.csv')

retainedIntron_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_retainedIntron_Acceptor_intron.csv')
retainedIntron_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_retainedIntron_Donor_intron.csv')
cassette_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_cassetteExon_Acceptor_intron.csv')
cassette_donor=read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_cassetteExon_Donor_intron.csv')

### Generating histograms ###
# altFive - Donor
ggplot(altFive_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 5' Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# AS - Acceptor
ggplot(altFive_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 5' Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# altThree - Donor
ggplot(altThree_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 3' Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# altThree - Acceptor
ggplot(altThree_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, fill="pink",color="pink") +
  labs(title="Stability distributions - Alternative 3' Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# Intron retention
ggplot(retainedIntron_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Retained intron Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
ggplot(retainedIntron_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.7, fill="pink",color="pink") +
  labs(title="Stability distributions - Retained intron Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")  +
  xlim(-100, 5)
# Constitutive - Donor
ggplot(CI_donor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence") +
  xlim(-100, 5)
# Constitutive - Acceptor
ggplot(CI_acceptor, aes(x=deltaG, fill=label)) +
  geom_histogram(position="identity", alpha=0.5, color="#619CFF", fill='#619CFF') +
  labs(title="Stability distributions - Constitutive Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

### Balanced data AS and CI ###
CI_balenced_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_constitutive_Donor.csv')
CI_balenced_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_constitutive_Acceptor.csv')

CI_balenced_donor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_cons_ALL_balanced_constitutive_Donor_intron.csv')
CI_balenced_acceptor = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_cons_ALL_balanced_constitutive_Acceptor_intron.csv')

# alternative Five
CI_altFive_constitutive = rbind(CI_balenced_acceptor, altFive_donor)
ggplot(CI_altFive_constitutive, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Acceptor vs Alt5 correct)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_alternative = rbind(CI_balenced_acceptor, altFive_acceptor)
ggplot(CI_altFive_alternative, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Acceptor vs Alt5 alternative)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_acceptor = rbind(CI_balenced_acceptor, altFive_acceptor)
ggplot(CI_altFive_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative 5 x Constitutive donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_acceptor = rbind(CI_balenced_acceptor, altFive_donor)
ggplot(CI_altFive_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative 5 - Donor splice site",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altFive_donor = rbind(CI_balenced_acceptor, altFive_acceptor)
ggplot(CI_altFive_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# alternative Three
CI_altThree_constitutive = rbind(CI_balenced_donor, altThree_acceptor)
ggplot(CI_altThree_constitutive, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Constitutive (CS Donor vs Alt3 correct)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_alternative = rbind(CI_balenced_donor, altThree_donor)
ggplot(CI_altThree_alternative, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative (CS Donor vs Alt3 alternative)",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_donor = rbind(CI_balenced_acceptor, altThree_donor)
ggplot(CI_altThree_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Alternative 3' - Acceptor Splice site",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_acceptor = rbind(CI_balenced_donor, altThree_acceptor)
ggplot(CI_altThree_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - CS donor x AltThree acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_altThree_acceptor = rbind(CI_balenced_donor, altThree_donor)
ggplot(CI_altThree_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - CS donor x AltThree donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# IR
CI_IR_acceptor = rbind(retainedIntron_acceptor, CI_balenced_acceptor)
ggplot(CI_IR_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_IR_donor = rbind(CI_balenced_donor, retainedIntron_donor)
ggplot(CI_IR_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Intron retention Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# Cassette
cassette_acceptor_balanced = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_cassetteExon_Acceptor_intron.csv')
cassette_donor_balanced = read.csv('/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/hg18_ALL_balanced_cassetteExon_Donor_intron.csv')
CI_cassette_aceptor = rbind(CI_balenced_acceptor, cassette_donor_balanced)
ggplot(CI_cassette_aceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exons - Acceptor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

CI_cassette_donor = rbind(CI_balenced_donor, cassette_acceptor_balanced)
ggplot(CI_cassette_donor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(title="Stability distributions - Cassette exon - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

### All together ### menos o cssette
# Acceptor
all_together_acceptor= rbind(CI_balenced_acceptor, altThree_acceptor, altFive_acceptor, retainedIntron_acceptor )
ggplot(all_together_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.3) +
  labs(title="Stability distributions - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")

# Donor
all_together_donor= rbind(CI_balenced_donor, altThree_donor, altFive_donor, retainedIntron_donor)
ggplot(all_together_acceptor, aes(x=deltaG, color=label,fill=label)) +
  geom_histogram(position="identity", alpha=0.3) +
  labs(title="Stability distributions - Donor",x="Stability of RNA secondary structure (kcals/mol)", y = "Occurence")



