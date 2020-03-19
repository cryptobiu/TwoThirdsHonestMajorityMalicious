#! /bin/bash
for i in $(seq "${1}" 1 "${2}");
do	
	./TwoThirdsHonestMajorityMalicious -partyID "${i}" -partiesNumber "${3}" -inputFile "${4}" -outputFile output.txt \
	 -circuitFile "${5}" -fieldType "${6}" -partiesFile "${7}" -isHonest "${8}"  -internalIterationsNumber "${9}" \
	   -isPRF "${10}" -numberOfTriples "${12}" &
	echo "Running $i..."
done
