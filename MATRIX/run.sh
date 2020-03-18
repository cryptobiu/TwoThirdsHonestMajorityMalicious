#!/usr/bin/env bash

./TwoThirdsHonestMajorityMalicious -partyID "${1}" -partiesNumber "${2}" -inputFile "${3}" -outputFile output.txt \
-circuitFile "${4}" -fieldType "${5}" -partiesFile "${6}" -isHonest "${7}" -internalIterationsNumber "${8}" \
-isPRF "${9}" -numberOfTriples "${10}"