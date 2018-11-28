#!/bin/bash

SERVERPATH="https://download.geofabrik.de/asia"
SUFFIX="-latest-free.shp.zip"

for COUNTRY in vietnam cambodia thailand india
do
    FILENAME=${COUNTRY}${SUFFIX}
    wget -N ${SERVERPATH}/${FILENAME}
    mkdir -p ${COUNTRY}
    unzip -u -d ${COUNTRY} ${FILENAME}
done
