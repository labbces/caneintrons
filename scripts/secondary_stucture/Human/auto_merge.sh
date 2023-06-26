#!/bin/bash

LISTA=("altFinish" "altFivePrime" "altPromoter" "altThreePrime" "atacIntron" "bleedingExon" "cassetteExon"  "retainedIntron" "strangeSplice")
CAMINHO=/data/AG-BioInf/guestwalther3/human_genome

for AS_TYPE in ${LISTA[@]}
do
        echo ${AS_TYPE}
	cat  hg18_chr*_${AS_TYPE}_Acceptor_intron.fa >>  hg18_ALL_${AS_TYPE}_Acceptor_intron.fa
done