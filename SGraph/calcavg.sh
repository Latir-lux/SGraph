#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Input a valid path."
    exit
fi

directory="$1"
filelist=`ls ${directory}`
rowcount=0
filecount=0

for file in $filelist
do
    if [[ $file =~ \.txt$ ]]; then
        row=`awk 'END{print NR}' ${directory}${file}`
        filecount=`expr ${filecount} + 1`
        rowcount=`expr ${rowcount} + ${row}`
    fi
done

echo `expr ${rowcount} / ${filecount}`
