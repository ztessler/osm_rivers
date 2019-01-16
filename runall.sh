#!/bin/bash

DELTAs=(Mekong Chao_Phraya Godavari)
OSMrivers=(vietnam:cambodia thailand india)
#RESs=(30min 10min 06min 03min 01min)
RESs=(10min 06min 03min)

for RES in ${RESs[@]}
do
    for IDX in ${!DELTAs[@]}
    do
        export DELTA=${DELTAs[$IDX]}
        export STNres=$RES
        export OSMriver=${OSMrivers[$IDX]}
        scons
    done

    # build SSEA bifur file
    mkdir -p output/SSEA/${RES}
    find output/ -path output/SSEA -prune -o -name "*_${RES}_bifurcations.csv" -execdir cat '{}' \+ > output/SSEA/${RES}/SSEA_${RES}_bifurcations.csv
done



