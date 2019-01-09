#!/bin/bash

DELTAs=(Mekong Chao_Phraya Godavari)
OSMrivers=(vietnam:cambodia thailand india)
#RESs=(30min 10min 06min 03min 01min)
RESs=(03min 10min 06min)

for IDX in ${!DELTAs[@]}
do
    for RES in ${RESs[@]}
    do
        export DELTA=${DELTAs[$IDX]}
        export STNres=$RES
        export OSMriver=${OSMrivers[$IDX]}

        # first for each delta as seperate domain
        export DOMAIN=${DELTAs[$IDX]}
        scons

        # then for each delta in SSEA domain
        export DOMAIN="SSEA"
        scons
    done
done


