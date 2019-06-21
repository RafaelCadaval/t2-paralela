#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>
#include <omp.h>
#include <stdbool.h>

void printMatrix(int size, int matrix[size][SIZE]) {
    int row, columns;
    for (row=0; row<size; row++)
    {
        for(columns=0; columns<SIZE; columns++)
            {
            printf("%d     ", matrix[row][columns]);
            }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    bool MASTER;
    int my_rank;       // Identificador deste processo
    const int MASTER_RANK = 0;
    int proc_n;        // Numero de processos disparados pelo usuario na linha de comando (np)   
    int hostsize;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    char master_hostname[MPI_MAX_PROCESSOR_NAME];    
    bool application_finished = false;
    double elapsed_time;

    MPI_Status status;

    // ************* TAGS *************
    // Possible tags for the master to process
    const int MASTER_HOSTNAME_TAG = 0;
    const int MATRIX_B_BCAST_TAG = 1;
    const int LINES_TO_PROCESS_TAG = 2;
    const int STOP_WORKER = 3;
    // ...

    MPI_Init(&argc , &argv); // funcao que inicializa o MPI, todo o codigo paralelo esta abaixo
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // pega o numero do processo atual (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &proc_n);  // pega informacao do numero de processos (quantidade total)
    MPI_Get_processor_name(hostname, &hostsize);

    // Checks if the number of process is a valid one (>1)
    if(proc_n == 1) {
        printf("Número de processos inválido. Finalizando execução.\n");
        MPI_Finalize();
        exit(0);
    }

    if (my_rank == MASTER_RANK)
        MASTER = true;
    else
        MASTER = false;

    /*
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
    */

    // Starts setup timer
    // elapsed_time = - MPI_Wtime();

    if(MASTER) {    
        // *** Master ***

        MPI_Get_processor_name(master_hostname, &hostsize);
        // Broadcasting the master's hostname to all workers
        MPI_Bcast(&master_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MASTER_HOSTNAME_TAG, MPI_COMM_WORLD);

        // Inicialização das matrizes
        int matrixA[SIZE][SIZE], matrixB[SIZE][SIZE];//, matrixRes[3][3];
        
        // Matriz A
        int k=1;
        int i, j;
        for (i=0 ; i<SIZE; i++) {
            for (j=0 ; j<SIZE; j++) {
                if (k%2==0)
                    matrixA[i][j] = -k;
                else
                    matrixA[i][j] = k;
            }
            k++;
        }    

        // Matriz B
        k=1; i=0; j=0;
        for (i=0 ; i<SIZE; i++) {
            for (j=0 ; j<SIZE; j++) {
                if (k%2==0)
                    matrixB[i][j] = -k;
                else
                    matrixB[i][j] = k;
                }
            k++;
        }

        printf("*** MASTER MATRIXES ***\n");
        printf("***    A MATRIX    ***\n");
        printMatrix(SIZE, matrixA);
        printf("**********************\n");
        printf("**********************\n");
        printf("***    B MATRIX    ***\n");
        printMatrix(SIZE, matrixB);
        printf("**********************\n");
        printf("**********************\n\n");

        // Broadcast matrix B to all workers
        MPI_Bcast(&matrixB, SIZE*SIZE, MPI_INT, MATRIX_B_BCAST_TAG, MPI_COMM_WORLD);

        // while(!application_finished) {
        //     // MPI_Recv(buffer, count, data_type, source, tag, comm, status)
        //     // MPI_Send(buffer, count, data_type, dest, tag, comm)

        //     MPI_Recv();
        // }
        
    } else {
        // *** Worker ***

        // Waiting for the master's hostname broadcast
        MPI_Bcast(&master_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MASTER_HOSTNAME_TAG, MPI_COMM_WORLD);
        // Setting up the number of threads (with two options):
        //  - num_threads = number of the node's processors, if the master isn't running on the same node;
        //      OR
        //  - num_threads = number of the node's processors - 1, if the master is running on the same node.
        int num_threads = omp_get_num_procs();
        if(strcmp(hostname, master_hostname)==0) {
            --num_threads;
        }
        omp_set_num_threads(num_threads);

        printf("\nNumber of threads: %d\n", num_threads);
        
        int matrixBAux[SIZE][SIZE];
        MPI_Bcast(&matrixBAux, SIZE*SIZE, MPI_INT, MATRIX_B_BCAST_TAG, MPI_COMM_WORLD);

        printf("** Bcast matrix before pFor **\n");
        printMatrix(SIZE, matrixBAux);
        printf("**********************\n");

        int i;
        #pragma omp parallel for
        for(i=0;i<num_threads;i++) {
            printf("** Worker %d matrix **\n", i);
            printMatrix(SIZE, matrixBAux);
            printf("**********************\n", i);
        }

        // MPI_Recv();

    }

    // Stops execution timer
    // elapsed_time += MPI_Wtime();
    // printf("Setup total time: %f seconds\n", elapsed_time);

    // Syncs the execution's start
    // MPI_Barrier(MPI_COMM_WORLD);
    // Starts execution timer
    // elapsed_time = - MPI_Wtime();

    // code

    // Stops execution timer
    // elapsed_time += MPI_Wtime();
    // printf("Execution total time: %f seconds\n", elapsed_time);

    MPI_Finalize();
    return 0;
}
