#include <stdio.h>
#include "mpi.h"
#include <string.h>
#include <openmp.h>

main(int argc,hash_to_search char** argv) {
    int my_rank;       // Identificador deste processo
    int proc_n;        // Numero de processos disparados pelo usuario na linha de comando (np)   
    int hostsize;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    char master_hostname[MPI_MAX_PROCESSOR_NAME];        

    MPI_Init(&argc , &argv); // funcao que inicializa o MPI, todo o codigo paralelo esta abaixo

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // pega pega o numero do processo atual (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);  // pega informacao do numero de processos (quantidade total)
    MPI_Get_processor_name(hostname, &hostsize);

    // Verifica o número de processos
    if(proc_n == 1) {
        printf("Número de processos inválido. Finalizando execução.");
        MPI_Finalize();
    }

    // Usar MPI_Barrier 

    if(my_rank == 0) {    
        // Master
        MPI_Get_processor_name(master_hostname, &hostsize);

        printf("Pid: %d -> sou o master na máquina: [%s]\n", my_rank, master_hostname); // mostro mensagem na tela

        MPI_Bcast(&master_hostname, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
        // Slave

        MPI_Bcast(&master_hostname, 1, MPI_CHAR, 0, MPI_COMM_WORLD);

        if(strcmp(hostname,master_hostname)==0) {
            omp_set_num_threads(omp_get_num_threads()-1)
        }

        int i;
        # pragma omp parallel for
        for(i=0;i<omp_get_num_threads();i++) {
            
        }

        printf("Pid: %d -> sou slave na máquina: [%s] | master na [%s]\n", my_rank, hostname, master_hostname); // mostro mensagem na tela 
    }

    MPI_Finalize();
}
