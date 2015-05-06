#!/bin/bash

run=$1

for file in *.list; do
    cat $file >> event.list
done

cp wn.xml wn.log

time root -q -b "rawmerge.C(\"wn.xml\",\"event.list\",\"filtered.root\")" &>filter.log
