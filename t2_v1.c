#include <stdio.h>
#include "mpi.h"

main(int argc, char** argv) {
    int my_rank;       // Identificador deste processo
    int proc_n;        // Numero de processos disparados pelo usuario na linha de comando (np)   
    int hostsize;
    char hostname[MPI_MAX_PROCESSOR_NAME];          

    MPI_Init(&argc , &argv); // funcao que inicializa o MPI, todo o codigo paralelo esta abaixo

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // pega pega o numero do processo atual (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);  // pega informacao do numero de processos (quantidade total)
    MPI_Get_processor_name(hostname, &hostsize);

    if(my_rank == 0) {
        // Master
        printf("Pid: %d -> sou o master na máquina: [%s]\n", my_rank, hostname); // mostro mensagem na tela 
    } else {
        // Slave
        printf("Pid: %d -> sou slave na máquina: [%s]\n", my_rank, hostname); // mostro mensagem na tela 
    }

    MPI_Finalize();
}
