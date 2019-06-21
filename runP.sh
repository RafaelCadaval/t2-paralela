git pull &&
mpicc parallel.c -o parallel -DSIZE=5 -std=c99 -fopenmp &&
mpirun -np 1 parallel
