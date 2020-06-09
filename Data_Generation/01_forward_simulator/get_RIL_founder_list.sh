#!/usr/bin/env bash

for file in RIL*.haps; do
    population=${file%.haps}
    cat $file | cut -f 7 | sort -u | grep -v "lineID" > ${population}.founders
done
