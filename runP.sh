mpicc parallel.c -o parallel -DSIZE=2000 -std=c99 -fopenmp &&
mpirun -np 2 parallel
