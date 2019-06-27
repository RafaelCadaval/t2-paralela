#!/bin/bash

i=100
while [ "$i" -le 2000 ]; do
    j=1
    while [ "$j" -le 10 ]; do
        mpicc parallel.c -o parallel -DSIZE=2000 -std=c99 -fopenmp &&
        echo "$i - $j " && 
        qsub batchjob-2/Batchjob-2-proc >> output.txt 
        j=$(( j + 1 ))
    done
    i=$(( i + 100 ))
done 