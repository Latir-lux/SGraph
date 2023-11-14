#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Input 2 valid path."
    exit
fi

directory1="$1"
directory2="$2"
filelist=`ls ${directory1}`

for file in $filelist
do
    if [[ $file =~ \.txt$ ]]; then
        `sort ${directory1}${file} | uniq > ${directory2}/uniq_${file}`
    fi
done
