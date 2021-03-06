#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>
#include <omp.h>
#include <stdbool.h>

#define MASTER_RANK 0

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
const int OFFSET_TAG = 5;
const int MATRIX_MULTIPLICATION_RESULT_TAG = 6;
const int HAS_ROWS_LEFT_TAG = 7;
const int CONTINUE_PROCESSING_TAG = 8;
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

    k=1;
    for (j=0; j<SIZE; j++) {
        for (i=0; i<SIZE; i++) {
            matrixRes[i][j] = 0;
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

void insertResMatrix(int start, int lines, int matrix[lines][SIZE]) {
    int i, j;
    for(i = 0; i < lines; i++) {
        for(j = 0; j < SIZE; j++) {
            matrixRes[i + start][j] = matrix[i][j];
        }
    }
}


void multiplyMatrix(int lines, int matrixA[lines][SIZE], int matrixB[SIZE][SIZE], int mres[lines][SIZE]) {
    int j, i, k;
    #pragma omp parallel for private(j, i, k)
    for (i=0 ; i<lines; i++) {
        for (j=0 ; j<SIZE; j++) {
            mres[i][j] = 0;
            for (k=0 ; k<SIZE; k++) {
                mres[i][j] += matrixA[i][k] * matrixB[k][j];
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

int main(int argc, char** argv) {
    // *** VARIABLES ***
    bool IS_MASTER;
    int my_rank; // Process ID
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

        num_threads = omp_get_num_procs();
        if(strcmp(hostname, master_hostname)==0) {
            --num_threads;
        }
        omp_set_num_threads(num_threads);
        MAX_NUMBER_OF_LINES = num_threads;

        MPI_Bcast(&matrixBAux, SIZE*SIZE, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
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

        int offset = 0;
        int rows_remaining = SIZE;
        bool should_continue = true;
        while(should_continue) {
            MPI_Status status;
            int throwaway_buffer;
            
            if(rows_remaining == 0) {
                int i;
                for(i=1;i<num_proc;i++) {
                    MPI_Recv(&throwaway_buffer, 1, MPI_INT, MPI_ANY_SOURCE, CONTINUE_PROCESSING_TAG, MPI_COMM_WORLD, &status);
                    throwaway_buffer = 0;
                    MPI_Send(&throwaway_buffer, 1, MPI_INT, status.MPI_SOURCE, CONTINUE_PROCESSING_TAG, MPI_COMM_WORLD);
                }
                should_continue = false;
                continue;
            } else { 
                MPI_Recv(&throwaway_buffer, 1, MPI_INT, MPI_ANY_SOURCE, CONTINUE_PROCESSING_TAG, MPI_COMM_WORLD, &status);
                throwaway_buffer = 1;
                MPI_Send(&throwaway_buffer, 1, MPI_INT, status.MPI_SOURCE, CONTINUE_PROCESSING_TAG, MPI_COMM_WORLD);
                int num_rows;
                MPI_Recv(&num_rows, 1, MPI_INT, MPI_ANY_SOURCE, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, &status);

                if(rows_remaining < num_rows) {
                    num_rows = rows_remaining;
                } 

                int test_matrix[num_rows][SIZE];
                copyMatrix(offset, num_rows, matrixA, test_matrix);
                
                MPI_Send(&offset, 1, MPI_INT, status.MPI_SOURCE, OFFSET_TAG, MPI_COMM_WORLD);
                MPI_Send(&num_rows, 1, MPI_INT, status.MPI_SOURCE, NUMBER_OF_ROWS_TAG, MPI_COMM_WORLD);
                MPI_Send(&test_matrix, num_rows * SIZE, MPI_INT, status.MPI_SOURCE, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);
                int matrixResCut[num_rows][SIZE];
                MPI_Recv(&matrixResCut, num_rows * SIZE, MPI_INT, MPI_ANY_SOURCE, MATRIX_MULTIPLICATION_RESULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                insertResMatrix(offset, num_rows, matrixResCut);
            
                rows_remaining -= num_rows;
                offset += num_rows;
            }
        }
        
        if(verifyResult()) {
            // verifyResult will return 0 if the multiplication is correct and 1 if it isn't
            printf("Matrix multiplication is not correct. Terminating execution.\n");
            MPI_Finalize();
            exit(0);
        } 

        // Stops execution timer
        execution_elapsed_time += MPI_Wtime();
        printf("\n\n*************************************\n");
        printf("Execution total time: %f seconds\n", execution_elapsed_time);
        printf("*************************************\n\n");

        MPI_Finalize();
    } else {
        // *** Worker process ***
        
        int has_rows_left = 1;
        while(has_rows_left) {
            int throwaway_buffer = 0;
            MPI_Send(&throwaway_buffer, 1, MPI_INT, MASTER_RANK, CONTINUE_PROCESSING_TAG, MPI_COMM_WORLD);
            MPI_Recv(&throwaway_buffer, 1, MPI_INT, MASTER_RANK, CONTINUE_PROCESSING_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            if(!throwaway_buffer) {  
                has_rows_left = throwaway_buffer;
                continue; 
            }

            MPI_Send(&MAX_NUMBER_OF_LINES, 1, MPI_INT, MASTER_RANK, NEED_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD);

            int offset;
            MPI_Recv(&offset, 1, MPI_INT, MASTER_RANK, OFFSET_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int num_rows;
            MPI_Recv(&num_rows, 1, MPI_INT, MASTER_RANK, NUMBER_OF_ROWS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int recv_matrix[num_rows][SIZE];
            MPI_Recv(&recv_matrix, num_rows * SIZE, MPI_INT, MASTER_RANK, SENDING_LINES_TO_PROCESS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int res[num_rows][SIZE];

            multiplyMatrix(num_rows, recv_matrix, matrixBAux, res);

            MPI_Send(&res, num_rows * SIZE, MPI_INT, MASTER_RANK, MATRIX_MULTIPLICATION_RESULT_TAG, MPI_COMM_WORLD);
        }
        MPI_Finalize();
    }

    return 0;
}
