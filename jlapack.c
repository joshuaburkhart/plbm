#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <clapack.h>

#define L 4

void output(double *matrix,int m,int n);
double* eye(int n);
double* kron(double *A,int ma,int na,double *B,int mb,int nb);
double* tran(double *A,int ma,int na);
double* 

int main(void){

  double *result;

  printf("---------------identity\n");

  result = eye(L);
  output(result,L,L);

  printf("---------------kronecker\n");

  double A[4]={1,2,3,4};
  double B[4]={0,5,6,7};
  result=kron(A,2,2,B,2,2);
  output(result,4,4);

  printf("---------------transpose\n");

  double C[6]={1,2,3,4,5,6};
  result=tran(C,2,3);
  output(result,3,2);
 
  free(result);
  return 0;
}

double* eye(int n){

  double *iden;
  iden = (double *) malloc(n * n * sizeof(double));
  memset(iden,0,n*n*sizeof(double));
  int i;
  int j;
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i==j){
        *(iden+(i*n+j))=1;
      }
    }
  }
  return iden;
}

double* kron(double *A,int ma,int na,double *B,int mb,int nb){

  double* k;
  k = (double *) malloc(ma*mb*na*nb*sizeof(double));
  int i;
  int j;
  for(i=0;i<ma*mb;i++){
    for(j=0;j<na*nb;j++){
      double val=*(A+((i/ma)*na+j/na)) * *(B+((i%mb)*nb+j%nb));
      *(k+((i*na*nb)+j))=val;
    }
  }
  return k;
}

double* tran(double *A,int ma,int na){

  double* t;
  t = (double *) malloc(ma*na*sizeof(double));
  int i;
  int j;
  for(i=0;i<ma;i++){
    for(j=0;j<na;j++){
      *(t+(i+j*ma))=*(A+(i*na+j));
    }
  }
  return t;
}

void output(double *matrix,int m,int n){

  int i;
  int j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
        printf("%f ",*(matrix+(i*n+j)));
    }   
    printf("\n");
  }
}
