#!/bin/bash

RESs=(10min 06min 03min 01min)

for RES in ${RESs[@]}
do
        export STNres=$RES
        scons -j 4

    # build SSEA bifur file
    mkdir -p output/SSEA/${RES}
    find output/ -path output/SSEA -prune -o -name "*_${RES}_bifurcations.csv" -execdir cat '{}' \+ > output/SSEA/${RES}/SSEA_${RES}_bifurcations.csv
done



