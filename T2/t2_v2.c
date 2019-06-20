#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>
#include <omp.h>

main(int argc, char** argv) {
    int my_rank;       // Identificador deste processo
    int proc_n;        // Numero de processos disparados pelo usuario na linha de comando (np)   
    int hostsize;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    char master_hostname[MPI_MAX_PROCESSOR_NAME];        

    MPI_Init(&argc , &argv); // funcao que inicializa o MPI, todo o codigo paralelo esta abaixo

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // pega pega o numero do processo atual (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);  // pega informacao do numero de processos (quantidade total)
    MPI_Get_processor_name(hostname, &hostsize);

 /*   // Verifica o número de processos
    if(proc_n == 1) {
        printf("Número de processos inválido. Finalizando execução.\n");
        MPI_Finalize();
	exit(0);
    }
*/
    // Inicialização das matrizes
    int matrixA[3][3], matrixB[3][3];
//, matrixRes[3][3];
    
    // Matriz A
    int k=1;\
    int i, j;
    for (i=0 ; i<3; i++) {
        for (j=0 ; j<3; j++) {
	  if (k%2==0)
               matrixA[i][j] = -k;
	  else
               matrixA[i][j] = k;
        }
        k++;
    }    

    // Matriz B
    k=1;
    i=0; 
    j=0;
    for (i=0 ; i<3; i++) {
        for (j=0 ; j<3; j++) {
          if (k%2==0)
               matrixB[i][j] = -k;
          else
               matrixB[i][j] = k;
        }
        k++;
    }

    printf("MATRIZES:\n");
    printf("[A]\n");
    int row, columns;
    for (row=0; row<3; row++) {
        for(columns=0; columns<3; columns++)
            {
                printf("%d     ", matrixA[row][columns]);
            }
        printf("\n");
    }
    printf("\n-------------\n");
    printf("[B]\n");
    for (row=0; row<3; row++) {
        for(columns=0; columns<3; columns++)
            {
                printf("%d     ", matrixB[row][columns]);
            }
        printf("\n");
    }
    printf("\n-------------\n");


    // Usar MPI_Barrier 
/*
    if(my_rank == 0) {    
        // Master
        MPI_Get_processor_name(master_hostname, &hostsize);

        printf("Master running on: [%s]\n", master_hostname);

        MPI_Bcast(&master_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
        // Worker
        MPI_Bcast(&master_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);
	int num_threads = omp_get_num_procs();
        if(strcmp(hostname, master_hostname)==0) {
            --num_threads;
      	}

	omp_set_num_threads(num_threads);
        int i;
        #pragma omp parallel for
        for(i=0;i<omp_get_num_threads();i++) {
            printf("Worker #%d running on [%s] | master on [%s]\n", i, hostname, master_hostname);
        } 
    }
*/
    MPI_Finalize();
}


