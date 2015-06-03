/* Prof.: Silvana Rossetto */
/* Introdução ao MPI (biblioteca de troca de mensagens) */
/* Mostra como usar as funções de distribuição e coleta de dados distribuídos */

/* Compile:  mpicc -O -o produtoMatrizVetor produtoMatrizVetor.c */
/* Run:      mpirun -np <number of processes> ./produtoMatrizVetor <tamanho do vetor> <debug_flag>*/
/* Exm:      mpirun -np 2 ./produtoMatrizVetor 1000*/
/* --------------------------------------------------------------------------------------------*/

#include "mpi.h"
#include "timer.h"
#include <stdio.h>
#include <stdlib.h>

typedef int bool;
#define true 1
#define false 0

#define SEED 1
#define DEBUG 1

int** InitMatriz(int n, bool cheia);
int** InitMatrizCheia(int tamanhoMatriz);
int** InitMatrizVazia(int tamanhoMatriz);
//void PrintMatriz(int** matrint** InitMatriz(int n,bool cheia)iz,int tamanhoMatriz);
void PrintMatriz(int* matriz,int tamanhoMatriz);
void PrintVetor(int* vetor, int tamanhoVetor);

int* ProdutoSequencial(int* A, int* X, int tamanho){
    int* Y;
    Y=malloc(tamanho*sizeof(int));
    int i,j;
    for(i=0;i<tamanho;i++){
        Y[i]=0;
        for(j=0;j<tamanho;j++){
            Y[i]+=A[j+i*tamanho]*X[j];
        }
    }
    return Y;
}

bool ComparaVetores(int* A, int* B, int tamanho){
    int i;
    for(i=0;i<tamanho;i++){
        if(A[i]!=B[i]) return false;
    }
    return true;
}

int main(int argc, char *argv[]) {
   int rank, nprocs, i, j, ret, local_flag=0, total_flag=0;
   int *local_a=NULL, size, local_n;
   int *A=NULL,*X=NULL,*Y=NULL,*R=NULL;
   double startTimeParallel,endTimeParallel, elapsedTimeParallel;
   double startTimeSequential,endTimeSequential, elapsedTimeSequential;

   if(SEED) srand(SEED);
   else srand(time(NULL));

   //inicializa o MPI
   ret = MPI_Init(&argc, &argv);
   if (ret == MPI_SUCCESS){
     MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   } else return 0;

   //recebe o numero de elementos do vetor
   if(argc<2) {
      if(rank==0) printf("Digite: %s <numero de elementos do vetor>\n", argv[0]);
      MPI_Finalize();
      return 0;
   }
   size = atoi(argv[1]);

    X = malloc( size*sizeof( int ));
    if(X==NULL) local_flag = 1; //sinaliza que houve erro
   //aloca espaco para o vetor de entrada e de saida
   if (rank == 0) {
       A = malloc( size*size*sizeof( int ));
#ifdef DEBUG
        printf("A\n");
#endif
       if(A==NULL) local_flag = 1; //sinaliza que houve erro
       for(i=0;i<size*size;i++) A[i]=i;

       Y = malloc( size*sizeof( int ));
       if(Y==NULL) local_flag = 1; //sinaliza que houve erro
#ifdef DEBUG
        printf("Y:\n");
#endif

       for(i=0;i<size;i++) X[i]=rand() % 10;
   }
#ifdef DEBUG
    if(rank==0){
        printf("A:\n");
        PrintMatriz(A,size);
        printf("X:\n");
        PrintVetor(X,size);
    }
#endif
   //aloca espaco para o bloco do vetor em cada processo
   local_n = (size/nprocs)*size;
   local_a = (int*) malloc(local_n * sizeof(int));
   if (local_a == NULL) local_flag = 1; //sinaliza que houve falha

   MPI_Bcast(X,size,MPI_INT,0,MPI_COMM_WORLD);
   //trata a possibilidade de falha na alocacao de memoria
   MPI_Allreduce (&local_flag, &total_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if(total_flag > 0) {
      printf("Erro de alocacao de memoria\n");
      MPI_Finalize(); //encerra todos os processos
      return 0;
   }
    #ifndef DEBUG
        if(rank==0) GET_TIME(startTimeParallel);
    #endif
   //divide o vetor em blocos eentrega um bloco por processo
   MPI_Scatter(A, local_n, MPI_INT, local_a, local_n, MPI_INT, 0, MPI_COMM_WORLD);
   int minhasLinhas = size/nprocs;
   int* acumulador = malloc(minhasLinhas*sizeof(int));
   //faz o processamento em cada processo
   int offset;
#ifdef DEBUG
printf("minhasLinhas=%d\n",minhasLinhas);
printf("[%d]:for(i =%d;i<%d;i++)\n",rank,0,minhasLinhas);
#endif
   for(i =0;i<minhasLinhas;i++){
       acumulador[i]=0;
       for (j = 0; j < size; j++){
          offset = j+i*size;
#ifdef DEBUG
printf("[%d]:acumulador[%d]+=local_a[%d]*X[%d];\n",rank,i,offset,j);
#endif
          acumulador[i]+=local_a[offset]*X[j];
       }
       #ifdef DEBUG
       printf("p:[%d] acumulador[%d]==%d\n",rank,i, acumulador[i]);
       #endif
   }
   if(rank==0){
#ifdef DEBUG
printf("%d*%d!=%d\n",minhasLinhas,nprocs,size);
#endif
    if(minhasLinhas*nprocs!=size){
        for(i=minhasLinhas*nprocs;i<size;i++){
           Y[i]=0;
           for (j = 0; j < size; j++){
              offset = j+i*size;
#ifdef DEBUG
printf("[%d]:Y[%d]+=A[%d]*X[%d];\n",rank,i,offset,j);
#endif
              Y[i]+=A[offset]*X[j];
           }
#ifdef DEBUG
if(DEBUG) printf("p:[%d] Y[%d]==%d\n",rank,i, Y[i]);
#endif // DEBUG
        }
    }
   }
      //MPI_Finalize(); //encerra todos os processos
      //return 0;
   //coleta os blocos processados por cada processo
   MPI_Gather(acumulador, minhasLinhas, MPI_INT, Y, minhasLinhas, MPI_INT, 0, MPI_COMM_WORLD);

    #ifndef DEBUG
        if(rank==0) GET_TIME(endTimeParallel);
    #endif
   //imprime o vetor resultante
   if (rank == 0) {
    #ifndef DEBUG
      GET_TIME(startTimeSequential);
      R=ProdutoSequencial(A,X,size);
      GET_TIME(endTimeSequential);
        elapsedTimeSequential = endTimeSequential - startTimeSequential;
        elapsedTimeParallel = endTimeParallel-startTimeParallel;
        printf("Sequential Time:%.6f\n",elapsedTimeSequential);
        printf("Parallel Time:%.6f\n",elapsedTimeParallel);
        printf("Speedup:%.6f\n",elapsedTimeSequential/elapsedTimeParallel);
    #endif // DEBUG
    #ifdef DEBUG
        printf("Y:\n");
        PrintVetor(Y,size);
        PrintVetor(R,size);
    #endif
      if(ComparaVetores(Y,R,size)) printf("Correto\n");
      else printf("Errado\n");

      free(A);
      free(Y);
   }
   free(X);
   free(local_a);
   free (acumulador);

   MPI_Finalize();
   return 0;
}

//void PrintMatriz(int** matriz,int tamanhoMatriz){
//    int i,j;
//    for(i=0;i<tamanhoMatriz;i++){
//        for(j=0;j<tamanhoMatriz;j++){
//            printf("%d ",matriz[i][j]);
//        }
//        printf("\n");
//    }
//}

void PrintMatriz(int* matriz,int tamanhoMatriz){
    int i,offset;
    if(tamanhoMatriz>10) return;
    for(i=0;i<tamanhoMatriz;i++){
        offset = tamanhoMatriz*i;
        PrintVetor(matriz+offset,tamanhoMatriz);
    }
}

void PrintVetor(int* vetor, int tamanhoVetor){
    int j;
    for(j=0;j<tamanhoVetor;j++){
        printf("%d ",vetor[j]);
    }
    printf("\n");
}

