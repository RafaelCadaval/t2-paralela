#!/bin/bash

i=100
while [ "$i" -le 2000 ]; do
    j=1
    while [ "$j" -le 10 ]; do
        echo "$i - $j " && qsub batchjob-2/Batchjob-2-proc >> output.txt 
        j=$(( j + 1 ))
    done
    i=$(( i + 100 ))
done 