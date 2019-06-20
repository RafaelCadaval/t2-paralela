git pull &&
mpicc parallel.c -o parallel -DSIZE=5 -std=c99 &&
mpirun -np 1 parallel
