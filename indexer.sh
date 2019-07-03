#!/bin/env bash

test

mkdir index
path=`pwd`
for file in `ls mapping`; do echo -e "$file\t$path/mapping\t$path/bac_analysis\t$path/sequences\t$path/primary_analysis" > index/$file; done
cat index/* > index.genome
