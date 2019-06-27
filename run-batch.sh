#!/bin/bash

for (( j=100; j<=2000; j+=100 ))
do
    for (( c=1; c<=10; c++ ))
    do
    mpicc parallel.c -o parallel -DSIZE=$j -std=c99 -fopenmp
    echo "$j - $c - ${qsub batchjob-2/Batchjob-2-proc}" >> output.txt
    done
done
