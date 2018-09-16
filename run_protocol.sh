#! /bin/bash
for i in `seq $1 1 $2`;
do	
	./MPCHonestMajorityNoTriplesYehuda -partyID $i -partiesNumber $3 -inputFile $4 -outputFile output.txt -circuitFile $5 -fieldType $6 -partiesFile $7 -isHonest 1  -internalIterationsNumber $8 &
	echo "Running $i..."
done
