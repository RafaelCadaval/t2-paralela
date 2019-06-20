git pull &&
mpicc sequencial.c -o sequencial -DSIZE=100 -std=c99 &&
mpirun -np 1 sequencial