#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Input a valid path."
    exit
fi

directory="$1"
filelist=`ls ${directory}`

for file in $filelist
do
    if [[ $file =~ \.txt$ ]]; then
        row=`awk 'END{print NR}' ${directory}${file}`
        if [ ${row} -gt 100000 ]; then
            echo ${row}
        else
            `rm -f ${directory}${file}`
        fi
    fi
done
