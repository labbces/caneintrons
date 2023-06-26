#!/bin/bash

for FILE in /home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/RNAfold_AS/*
do 
    FILENAME="$(basename -- $FILE)"
    tipo=$(echo "$FILENAME" | awk -F "_" '{ print $3 }')  
    FILENAME=${FILENAME/'RNAfold'/"csv"}
    FILENAME=/home/bia/sugarcane_introns_local/SecondaryStructure/Human_genome/csv/$FILENAME
    python3 -u "/home/bia/sugarcane_introns_local/caneintrons/scripts/secondary_stucture/RNAfold_analysis.py" -i ${FILE} -o ${FILENAME} -t $tipo
    echo "Output file $FILENAME generated"
done
