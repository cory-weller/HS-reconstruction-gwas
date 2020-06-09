#!/usr/bin/env bash

for file in *.haps; do
    population=${file%.haps}
    cat $file | cut -f 8 | sort -u | grep -v "lineID" > ${population}.founders
done
