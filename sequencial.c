// Arquivo: sequencial.c
// Autor    Roland Teodorowitsch
// Data:    28 ago. 2019

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// DADOS COMPARTILHADOS
int m1[SIZE][SIZE],m2[SIZE][SIZE],mres[SIZE][SIZE];
int l1, c1, l2, c2, lres, cres;

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

int initMatrixs() {
    int  i, j, k, id, p;
    l1 = c1 = SIZE;
    l2 = c2 = SIZE;
    if (c1 != l2) {
        fprintf(stderr,"Impossivel multiplicar matrizes: parametros invalidos.\n");
        return 1;
    }
    lres = l1;
    cres = c2;
    k=1;
    for (i=0 ; i<SIZE; i++) {
        for (j=0 ; j<SIZE; j++) {
        if (k%2==0)
                m1[i][j] = -k;
        else
                m1[i][j] = k;
        }
        k++;
    }
    k=1;
    for (j=0 ; j<SIZE; j++) {
        for (i=0 ; i<SIZE; i++) {
        if (k%2==0)
                m2[i][j] = -k;
        else
                m2[i][j] = k;
        }
        k++;
    }
    return 0;
}

int verifyResult() {
    int  i, j, k, id, p;
    for (i=0 ; i<SIZE; i++) {
      k = SIZE*(i+1);
      for (j=0 ; j<SIZE; j++) {
          int k_col = k*(j+1);
          if (i % 2 ==0) {
             if (j % 2 == 0) {
                if (mres[i][j]!=k_col)
                   return 1;
             }
             else {
                if (mres[i][j]!=-k_col)
                   return 1;
             }
          }
          else {
             if (j % 2 == 0) {
                if (mres[i][j]!=-k_col)
                   return 1;
             }
             else {
                if (mres[i][j]!=k_col)
                   return 1;
             }
          }
      } 
  }
  return 0;
}

void copyMatrix(int start, int lines, int m1[SIZE][SIZE], int matrix[lines][SIZE]) {
    int i, j;
    for(i = 0; i < lines; i++) {
        for(j = 0; j < SIZE; j++) {
            matrix[i] = m1[i + start][j];
        }
    }
}

void calculateMatrix(int lines, int m1[lines][SIZE], int matrix[lines][SIZE]) {
    int j, i, k;
    for (i=0 ; i<lines; i++) {
        for (j=0 ; j<SIZE; j++) {
            matrix[i][j] = 0;
            for (k=0 ; k<SIZE; k++) {
                matrix[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

int main(int argc, char *argv[]) {
  int    i, j, k, id, p, result;
  double elapsed_time;
  int matrix[5][SIZE];

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  if (id != 0) {
     MPI_Finalize();
     exit(0);
  }

  // INICIALIZA OS ARRAYS A SEREM MULTIPLICADOS
  result = initMatrixs();
  if(result != 0) {
    MPI_Finalize();
    return result;
  }

  // PREPARA PARA MEDIR TEMPO
  elapsed_time = - MPI_Wtime ();

  // REALIZA A MULTIPLICACAO
  for (i=0 ; i<lres; i++) {
      for (j=0 ; j<cres; j++) {
          mres[i][j] = 0;
          for (k=0 ; k<c1; k++) {
              mres[i][j] += m1[i][k] * m2[k][j];
          }
      }
  }

  // OBTEM O TEMPO
  elapsed_time += MPI_Wtime();

  // VERIFICA SE O RESULTADO DA MULTIPLICACAO ESTA CORRETO
  result = verifyResult();
  if(result != 0) {
      MPI_Finalize();
      return result;
  }

  printMatrix(SIZE, m1);
  printf("\n");
  printMatrix(SIZE, m2);
  printf("\n");
  printMatrix(SIZE, mres);
  printf("\n");

  int copyM1[5][SIZE];
  copyMatrix(5, 5, m1, copyM1);
  printMatrix(5, copyM1);

  printf("\n");
  printMatrix(SIZE, m1);
    

  // MOSTRA O TEMPO DE EXECUCAO
  printf("%lf",elapsed_time);
  MPI_Finalize();
  return 0;
}