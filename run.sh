git pull &&
mpicc sequencial.c -o sequencial -DSIZE=2000 -std=c99 &&
mpirun -np 1 sequencial