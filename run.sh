git pull;
mpicc sequencial.c -o sequencial -DSIZE=100;
mpirun -np 1 sequencial;