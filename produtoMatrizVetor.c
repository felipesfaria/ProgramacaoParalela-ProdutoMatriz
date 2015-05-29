/* Prof.: Silvana Rossetto */
/* Introdução ao MPI (biblioteca de troca de mensagens) */
/* Mostra como usar as funções de distribuição e coleta de dados distribuídos */

/* Compile:  mpicc -o produtoMatrizVetor produtoMatrizVetor.c */
/* Run:      mpirun -np <number of processes> ./produtoMatrizVetor */
/* Exm:      mpirun -np 2 ./produtoMatrizVetor 5 */
/* --------------------------------------------------------------------------------------------*/

#include "mpi.h"
#include "timer.h"
#include <stdio.h>
#include <stdlib.h>

typedef int bool;
#define true 1
#define false 0

#define SEED 1
//#define DEBUG 1

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
   int rank, nprocs, i, j, ret, local_flag=0, total_flag=0,DEBUG;
   int *a=NULL, *local_a=NULL, *b=NULL, size, local_n;
   int *A,*X,*Y,*R;
   double startTime,endTime, elapsedTime;

   if(SEED) srand(SEED);
   else srand(time(NULL));


   //inicializa o MPI
   ret = MPI_Init(&argc, &argv);
   if (ret == MPI_SUCCESS){
     MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   } else return 0;

   //recebe o numero de elementos do vetor
   if(argc<3) {
      if(rank==0) printf("Digite: %s <numero de elementos do vetor>\n", argv[0]);
      MPI_Finalize();
      return 0;
   }
   size = atoi(argv[1]);
   DEBUG = atoi(argv[2]);

    X = malloc( size*sizeof( int ));
   //aloca espaco para o vetor de entrada e de saida
   if (rank == 0) {
       A = malloc( size*size*sizeof( int ));
       for(i=0;i<size*size;i++) A[i]=rand() % 10;
       if(A==NULL) local_flag = 1; //sinaliza que houve erro
       Y = malloc( size*sizeof( int ));
       if(Y==NULL) local_flag = 1; //sinaliza que houve erro

       if(X==NULL) local_flag = 1; //sinaliza que houve erro
       for(i=0;i<size;i++) X[i]=rand() % 10;
   }
    if(rank==0){
        if(DEBUG){
            printf("A:\n");
            PrintMatriz(A,size);
            printf("X:\n");
            PrintVetor(X,size);
        }
    }
   //aloca espaco para o bloco do vetor em cada processo
   local_n = size/nprocs*size;
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
    if(rank==0&&!DEBUG) GET_TIME(startTime);
   //divide o vetor em blocos eentrega um bloco por processo
   MPI_Scatter(A, local_n, MPI_INT, local_a, local_n, MPI_INT, 0, MPI_COMM_WORLD);
   int minhasLinhas = size/nprocs;
   int* acumulador = malloc(minhasLinhas*sizeof(int));
   //faz o processamento em cada processo
   int offset;
   if(DEBUG) printf("minhasLinhas=%d\n",minhasLinhas);
   if(DEBUG) printf("[%d]:for(i =%d;i<%d;i++)\n",rank,0,minhasLinhas);
   for(i =0;i<minhasLinhas;i++){
       acumulador[i]=0;
       for (j = 0; j < size; j++){
          offset = j+i*size;
          if(DEBUG) printf("[%d]:acumulador[%d]+=local_a[%d]*X[%d];\n",rank,i,offset,j);
          acumulador[i]+=local_a[offset]*X[j];
       }
       if(DEBUG) printf("p:[%d] acumulador[%d]==%d\n",rank,i, acumulador[i]);
   }
   if(rank==0){
   if(DEBUG) printf("%d*%d!=%d\n",minhasLinhas,nprocs,size);
    if(minhasLinhas*nprocs!=size){
        for(i=minhasLinhas*nprocs;i<size;i++){
           Y[i]=0;
           for (j = 0; j < size; j++){
              offset = j+i*size;
              if(DEBUG) printf("[%d]:Y[%d]+=A[%d]*X[%d];\n",rank,i,offset,j);
              Y[i]+=A[offset]*X[j];
           }
           if(DEBUG) printf("p:[%d] Y[%d]==%d\n",rank,i, Y[i]);
        }
    }
   }
      //MPI_Finalize(); //encerra todos os processos
      //return 0;
   //coleta os blocos processados por cada processo
   MPI_Gather(acumulador, minhasLinhas, MPI_INT, Y, minhasLinhas, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank==0&&!DEBUG) GET_TIME(endTime);
   //imprime o vetor resultante
   if (rank == 0) {
      R=ProdutoSequencial(A,X,size);
    if(!DEBUG){
        elapsedTime = endTime-startTime;
        printf("Time:%.6f\n",elapsedTime);
    }
    if(DEBUG){
        printf("Y:\n");
        PrintVetor(Y,size);
        PrintVetor(R,size);
    }
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

int** InitMatrizCheia(int tamanhoMatriz){
    return InitMatriz(tamanhoMatriz,true);
}

int** InitMatrizVazia(int tamanhoMatriz){
    return InitMatriz(tamanhoMatriz,false);
}

int** InitMatriz(int n,bool cheia){
    int i,j;
    int** matriz;
    if (( matriz = malloc( n*sizeof( int* ))) == NULL )
    { /* error */ }
    for ( i = 0; i < n; i++ )
    {
        if (( matriz[i] = malloc( n*sizeof( int* ) )) == NULL )
        { /* error */ }
        if(cheia){
            for(j=0;j<n;j++){
                matriz[i][j]=(j+j+i);
            }
        }else{
            for(j=0;j<n;j++){
                matriz[i][j]=0;
            }
        }
    }
    return matriz;
}
