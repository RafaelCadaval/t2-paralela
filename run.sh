git pull &&
mpicc sequencial.c -o sequencial -DSIZE=5 -std=c99 &&
mpirun -np 1 sequencial