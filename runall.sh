#!/bin/bash

DELTAs=(Mekong Chao_Phraya Godavari)
OSMrivers=(vietnam:cambodia thailand india)
#RESs=(30min 10min 06min 03min 01min)
RESs=(10min 06min 03min)

for IDX in ${!DELTAs[@]}
do
    for RES in ${RESs[@]}
    do
        export DELTA=${DELTAs[$IDX]}
        export DOMAIN=${DELTAs[$IDX]}
        export STNres=$RES
        export OSMriver=${OSMrivers[$IDX]}
        scons
    done
done


