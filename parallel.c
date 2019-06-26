#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>
#include <omp.h>
#include <stdbool.h>

// SHARED DATA
int matrixA[SIZE][SIZE], matrixB[SIZE][SIZE], matrixBAux[SIZE][SIZE], matrixRes[SIZE][SIZE];
int l1, c1, l2, c2, lres, cres;
int num_proc;
int workers_finished = 0;
bool application_finished = false;

// ************* TAGS *************
// Possible tags for the master to process
const int MASTER_HOSTNAME_TAG = 0;
const int NEED_LINES_TO_PROCESS_TAG = 2;
const int SENDING_LINES_TO_PROCESS_TAG = 3;
const int NUMBER_OF_ROWS_TAG = 4;
const int STOP_WORKER_TAG = 10;
// ...






void printMatrix(int size, int matrix[size][SIZE]) {
    int row, columns;
    for (row=0; row<size; row++) {
        for (columns=0; columns<SIZE; columns++) {
            printf("%d\t", matrix[row][columns]);
        }
        printf("\n");
    }
}

int initMatrixes() {
    int i, j, k, id, p;
    // Verifying lines and columns sizes
    l1 = c1 = SIZE; // REFACTOR - can use only SIZE for everything
    l2 = c2 = SIZE;
    if (c1 != l2) {
        fprintf(stderr,"Impossible to multiply matrixes: invalid parameters.\n");
        MPI_Finalize();
        return 1;
    }

    lres = l1;
    cres = c2;
    k=1;
    for (i=0; i<SIZE; i++) {
        for (j=0; j<SIZE; j++) {
            if (k%2==0)
                matrixA[i][j] = -k;
            else
                matrixA[i][j] = k;
        }
        k++;
    }

    k=1;
    for (j=0; j<SIZE; j++) {
        for (i=0; i<SIZE; i++) {
            if (k%2==0)
                matrixB[i][j] = -k;
            else
                matrixB[i][j] = k;
        }
        k++;
    }
    return 0;
}

void copyMatrix(int start, int lines, int matrixA[SIZE][SIZE], int matrix[lines][SIZE]) {
    int i, j;
    for(i = 0; i < lines; i++) {
        for(j = 0; j < SIZE; j++) {
            matrix[i][j] = matrixA[i + start][j];
        }
    }
}

void multiplyMatrix(int lines, int matrixA[lines][SIZE], int matrix[lines][SIZE]) {
    int j, i, k;
    for (i=0; i<lines; i++) {
        for (j=0; j<SIZE; j++) {
            matrix[i][j] = 0;
            for (k=0; k<SIZE; k++) {
                matrix[i][j] += matrixA[i][k] * matrixBAux[k][j];
            }
        }
    }
}

int verifyResult() {
    int  i, j, k, id, p, k_col;
    for (i=0; i<SIZE; i++) {
        k = SIZE*(i+1);
        for (j=0; j<SIZE; j++) {
            k_col = k*(j+1);
            if (i % 2 ==0) {
                if (j % 2 == 0) {
                    if (matrixRes[i][j]!=k_col)
                        return 1;
                } else {
                    if (matrixRes[i][j]!=-k_col)
                        return 1;
                }
            } else {
                if (j % 2 == 0) {
                    if (matrixRes[i][j]!=-k_col)
                        return 1;
                } else {
                    if (matrixRes[i][j]!=k_col)
                        return 1;
                }
            }
        } 
    }
    return 0;
}

// void processStatus(MPI_Status status) {
//     // Extract status' tag
//     int tag = status.MPI_TAG;
//     int source = status.MPI_SOURCE;
//     int throwaway_buffer;
//     switch(tag) {
//     case NEED_LINES_TO_PROCESS_TAG:
//         // MPI_Recv(buffer, count, data_type, source, tag, comm, status)
//         MPI_Recv(throwaway_buffer, MAX_NUMBER_OF_LINES, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         // TEST CODE
//         int test_lines_buffer[1][SIZE];
//         copyMatrix(status.MPI_SOURCE, 1, matrixA, &test_lines_buffer);
//         // MPI_Send(&lines, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
//         MPI_Send(&test_lines_buffer, MAX_NUMBER_OF_LINES, MPI_INT, source, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
//         break;
//     case STOP_WORKER_TAG:
//         MPI_Recv(throwaway_buffer, MAX_NUMBER_OF_LINES, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         workers_finished++;
//         if(workers_finished == num_proc-1)
//             application_finished = true;
//         break;
//     default: 
//         printf("--> Couldn't parse request..\n");
//         break;
//     }
// }

int main(int argc, char** argv) {
    // *** VARIABLES ***
    bool IS_MASTER;
    int my_rank; // Process ID
    const int MASTER_RANK = 0;
    // int num_proc; // Number of process given by the user through the `np` clause
    int num_threads;
    int hostsize;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    char master_hostname[MPI_MAX_PROCESSOR_NAME];    
    double setup_elapsed_time, execution_elapsed_time;
    int MAX_NUMBER_OF_LINES;

    // MPI_Status status; // Additional information about the `receive` operation after it completes

    MPI_Init(&argc , &argv); // Initializes MPI; all parallel code is below this statement
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // Gets the current process number (rank)
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);  // Gets information about the number of process (total number)
    MPI_Get_processor_name(hostname, &hostsize); // Gets the node's name where the current process is running

    // Checks if the number of processes is a valid one (>1)
    if(!num_proc > 1) {
        printf("Number of process is invalid. Terminating execution.\n");
        MPI_Finalize();
        exit(0);
    }

    if (my_rank == MASTER_RANK)
        IS_MASTER = true;
    else
        IS_MASTER = false;

    if(IS_MASTER) {    
        // *** Master process ***

        // Starts setup timer
        setup_elapsed_time = - MPI_Wtime();

        MPI_Get_processor_name(master_hostname, &hostsize);
        // Broadcasting the master's hostname to all workers
        MPI_Bcast(&master_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MASTER_HOSTNAME_TAG, MPI_COMM_WORLD);

        // Initializing matrixes
        int result = initMatrixes();
        if(result != 0) {
            MPI_Finalize();
            return result;
        }

        // printf("*** MASTER MATRIXES ***\n");
        // printf("***    A MATRIX    ***\n");
        // printMatrix(SIZE, matrixA);
        // printf("**********************\n");
        // printf("**********************\n");
        // printf("***    B MATRIX    ***\n");
        // printMatrix(SIZE, matrixB);
        // printf("**********************\n");
        // printf("**********************\n\n");

        // Broadcast matrix B to all workers
        MPI_Bcast(&matrixB, SIZE*SIZE, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    } else {
        // *** Worker process ***

        // Waiting for the master's hostname broadcast
        MPI_Bcast(&master_hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);
        // Setting up the number of threads (with two options):
        //  - num_threads = number of the node's processors, if the master isn't running on the same node;
        //      OR
        //  - num_threads = number of the node's processors - 1, if the master is running on the same node.
        
        // printf("** Bcast m_hn **\n");
        // printf("%s\n", master_hostname);

        num_threads = omp_get_num_procs();
        if(strcmp(hostname, master_hostname)==0) {
            --num_threads;
        }
        omp_set_num_threads(num_threads);
        MAX_NUMBER_OF_LINES = num_threads;

        // printf("\nNumber of threads: %d\n", num_threads);

        MPI_Bcast(&matrixBAux, SIZE*SIZE, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

        // printf("** Bcast matrix before pFor **\n");
        // printMatrix(SIZE, matrixBAux);
        // printf("**********************\n");

        // int i;
        // #pragma omp parallel for
        // for(i=0;i<num_threads;i++) {
        //     printf("** Worker %d matrix **\n", i);
        //     // printMatrix(SIZE, matrixBAux);
        //     // printf("**********************\n", i);
        // }
    }

    // Syncs all process
    MPI_Barrier(MPI_COMM_WORLD);
    if(IS_MASTER) {
        // Stops setup timer
        setup_elapsed_time += MPI_Wtime();
        printf("\n\n*************************************\n");
        printf("Setup total time: %f seconds\n", setup_elapsed_time);
        printf("*************************************\n\n");
        // Starts execution timer
        execution_elapsed_time = - MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // ******************************************************************************************************
    // ******************************************* MAIN EXECUTION *******************************************
    // ******************************************************************************************************

    if(IS_MASTER) {
        // *** Master process ***
        
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        int buffer_count;
        MPI_Get_count(&status, MPI_INT, &buffer_count);
        if(status.MPI_TAG == NEED_LINES_TO_PROCESS_TAG) {
            int num_rows;
            MPI_Recv(&num_rows, buffer_count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("** Master received request.. **\n");
            int test_matrix[num_rows][SIZE];
            copyMatrix(status.MPI_SOURCE, num_rows, matrixA, test_matrix);
            printf("** Master sending matrix portion to Worker %d.. **\n", status.MPI_SOURCE);
            MPI_Send(&num_rows, num_rows, MPI_INT, status.MPI_SOURCE, NUMBER_OF_ROWS_TAG, MPI_COMM_WORLD);
            MPI_Send(&test_matrix, sizeof test_matrix, MPI_INT, status.MPI_SOURCE, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
        }







        /* GOING TO REFACTOR ALL THIS CODE (relax) */

        // MPI_Status status;
        // MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // int buffer_count;
        // MPI_Get_count(&status, MPI_INT, &buffer_count);
        // int throwaway_buffer = 0;

        // printf("--> 1º Buffer count: %d\n", buffer_count);

        // if(status.MPI_TAG == NEED_LINES_TO_PROCESS_TAG) {
        //     // MPI_Recv(buffer, count, data_type, source, tag, comm, status)
        //     MPI_Recv(&throwaway_buffer, buffer_count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     // TEST CODE
        //     int test_lines_buffer[1][SIZE];
        //     copyMatrix(status.MPI_SOURCE, 1, matrixA, test_lines_buffer);
        //     // MPI_Send(&lines, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
        //     MPI_Send(&test_lines_buffer, SIZE*SIZE, MPI_INT, status.MPI_SOURCE, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
        // } else if(status.MPI_TAG == STOP_WORKER_TAG) {
        //     MPI_Recv(&throwaway_buffer, buffer_count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //     workers_finished++;
        //     if(workers_finished == MAX_NUMBER_OF_LINES)
        //         application_finished = true;
        // } else {
        //     printf("--> Couldn't parse request..\n");
        // }



        /* HAVE TO FIX THIS LATER */

        // while(!application_finished) {
        //     // MPI_Recv(buffer, count, data_type, source, tag, comm, status)
        //     // MPI_Send(buffer, count, data_type, dest, tag, comm)
        //     // MPI_Probe(source, tag, comm, status);
        //     MPI_Status status;

        //     MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        //     int buffer_count;
        //     MPI_Get_count(&status, MPI_INT, &buffer_count);
        //     // processStatus(status, &lines);
        //     // processStatus(status);


        //     // Extract status' tag
        //     // int tag = status.MPI_TAG;
        //     int source = status.MPI_SOURCE;
        //     int throwaway_buffer = 0;

            // if(status.MPI_TAG == NEED_LINES_TO_PROCESS_TAG) {
            //     // MPI_Recv(buffer, count, data_type, source, tag, comm, status)
            //     MPI_Recv(&throwaway_buffer, buffer_count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //     // TEST CODE
            //     int test_lines_buffer[1][SIZE];
            //     copyMatrix(status.MPI_SOURCE, 1, matrixA, test_lines_buffer);
            //     // MPI_Send(&lines, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
            //     MPI_Send(test_lines_buffer, SIZE*SIZE, MPI_INT, status.MPI_SOURCE, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
            // } else if(status.MPI_TAG == STOP_WORKER_TAG) {
            //     MPI_Recv(&throwaway_buffer, buffer_count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //     workers_finished++;
            //     if(workers_finished == MAX_NUMBER_OF_LINES)
            //         application_finished = true;
            // } else {
            //     printf("--> Couldn't parse request..\n");
            // }



        // }
    } else {
        // *** Worker process ***

        // MPI_Send(&lines, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
        // MPI_Recv(&lines, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, status);

        printf("** Worker %d requesting lines **\n", my_rank);
        MPI_Send(&MAX_NUMBER_OF_LINES, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Probe(MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("** hello rank %d**\n", status.MPI_SOURCE);
        int buffer_count;
        MPI_Get_count(&status, MPI_INT, &buffer_count);
        if(status.MPI_TAG == NUMBER_OF_ROWS_TAG) {
            int num_rows;
            MPI_Recv(&num_rows, buffer_count, MPI_INT, MASTER_RANK, NUMBER_OF_ROWS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("** num_rows %d received**\n", num_rows);
            MPI_Probe(MASTER_RANK, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, &status);
            printf("** hello? **\n");
            MPI_Get_count(&status, MPI_INT, &buffer_count);
            int recv_matrix[num_rows][SIZE];
            MPI_Recv(&recv_matrix, buffer_count, MPI_INT, MASTER_RANK, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("** Worker %d received rows **\n", my_rank);
            printMatrix(num_rows, recv_matrix);
        }








        /* REFACTORING ALL THIS */

        // MPI_Status status;
        // printf("** Worker %d requesting lines **\n", my_rank);
        // int throwaway_buffer = 0;

        // printf("--> 2º Buffer count: %d\n", MAX_NUMBER_OF_LINES);
        // MPI_Send(&throwaway_buffer, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);

        
        // MPI_Probe(MASTER_RANK, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, &status);
        // int buffer_count;
        // MPI_Get_count(&status, MPI_INT, &buffer_count);
        // printf("--> 3º Buffer count: %d\n", buffer_count);
        // MPI_Recv(&lines, buffer_count, MPI_INT, MASTER_RANK, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // printf("** Lines given to Worker %d:\n", my_rank);
        // printMatrix(1, lines);
        // MPI_Send(&throwaway_buffer, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, STOP_WORKER_TAG, MPI_COMM_WORLD);





        /* HAVE TO FIX THIS LATER */

        // int i;
        // MPI_Status status;
        // #pragma omp parallel for private(lines, status)
        // for(i=0; i<MAX_NUMBER_OF_LINES; i++) { 
        //     printf("** Worker %d requesting lines **\n", i+1);
        //     int throwaway_buffer = 0;

        //     MPI_Send(&throwaway_buffer, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
            
        //     MPI_Probe(MASTER_RANK, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, &status);
        //     int buffer_count;
        //     MPI_Get_count(&status, MPI_INT, &buffer_count);
        //     MPI_Recv(&lines, buffer_count, MPI_INT, MASTER_RANK, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //     printf("** Lines given to Worker %d\n", i);
        //     printMatrix(1, lines);
        //     MPI_Send(&throwaway_buffer, MAX_NUMBER_OF_LINES, MPI_INT, MASTER_RANK, STOP_WORKER_TAG, MPI_COMM_WORLD);
        // }
    }

    // Stops execution timer
    // execution_elapsed_time += MPI_Wtime();
    // printf("\n\n*************************************\n");
    // printf("Execution total time: %f seconds\n", execution_elapsed_time);
    // printf("*************************************\n\n");

    MPI_Finalize();
    return 0;
}
