#!/bin/bash

SPP=human_genome
INPATH=/home/mpimp-golm.mpg.de/guestwalther3/data/${SPP}/RNAfold_results/*.RNAfold

for FILE in ${INPATH}
do
	# conda activate /home/mpimp-golm.mpg.de/guestwalther3/progs/cane_introns 
	FILENAME="$(basename -- $FILE)"
	FILENAME=${FILENAME/'.RNAfold'/'.csv'}
	echo "working on ${FILENAME}..."
	EVENT=$(echo $FILENAME | cut -d "_" -f 4)
	python3 /home/mpimp-golm.mpg.de/guestwalther3/caneintrons/scripts/secondary_stucture/RNAfold_analysis.py -i ${FILE} -o ${FILENAME} -t ${EVENT}
done
