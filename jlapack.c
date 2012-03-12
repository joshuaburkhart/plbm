#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <clapack.h>

#define L 4

void output(double *matrix,int m,int n);
double* eye(int n);
double* kron(double *A,int ma,int na,double *B,int mb,int nb);
double* tran(double *A,int m,int n);
double* ones(int m,int n);
double matrx_det(double *A,int n);
void dgesv_(const int *N, const int *nrhs,double *A,const int *lda,int *ipiv,double *b,const int *ldb,int *info);
void dgels_(const char *trans,const int *M,const int *N,const int *nrhs,double *A,const int *lda,double *b,const int *ldb,double *work,const int * lwork,int *info);
double* array_pow(double d,double *A,int m,int n);
double* array_mlt(double *A,int m,int n,double *B);
double* matrx_mlt(double d,double *A,int ma,int na);
double* matrx_mlt2(double *A,int ma,int na,double *B,int mb,int nb);
double* matrx_sub(double d,double *A,int ma,int na);
double* matrx_sub2(double *A,int ma,int na,double *B,int mb,int nb);
double* array_rdv(double *A,int ma,int na,double d);

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

  printf("---------------ones\n");

  result=ones(5,7);
  output(result,5,7);
 
  printf("---------------determinant\n");

  double d=matrx_det(A,2);
  printf("%f\n",d);

  printf("---------------array_pow\n");

  result=array_pow(5.00,A,2,2);
  output(result,2,2);

  printf("---------------array_mlt\n");

  result=array_mlt(A,2,2,B);
  output(result,2,2);

  printf("---------------matrx_mlt\n");

  result=matrx_mlt(5.00,A,2,2);
  output(result,2,2);

  printf("---------------matrx_mlt2\n");

  result=matrx_mlt2(A,2,2,B,2,2);
  output(result,2,2);

  free(result);
  return 0;
}

double* matrx_mlt2(double *A,int ma,int na,double *B,int mb,int nb){

  double *pdct;
  pdct=(double *) malloc(ma*nb*sizeof(double));
  int i;
  int j;
  int k;
  for(i=0;i<ma;i++){
    for(j=0;j<nb;j++){
      double sum=0;
      for(k=0;k<na;k++){
        sum+=*(A+(i*na+k)) * *(B+(j+k*mb));
      }
      *(pdct+(i*nb+j))=sum;
    }
  }
  return pdct; 
}

double* matrx_mlt(double d,double *A,int m,int n){

  double *pdct;
  pdct=(double *) malloc(m*n*sizeof(double));
  int i;
  int j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
        *(pdct+(i*n+j))=*(A+(i*n+j)) * d;
    }
  }
  return pdct;
}

double* array_mlt(double *A,int m,int n,double *B){

  double *pdct;
  pdct=(double *) malloc(m*n*sizeof(double));
  int i;
  int j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
        *(pdct+(i*n+j))=*(A+(i*n+j)) * *(B+(i*n+j));
    }
  }
  return pdct;
}

double* array_pow(double d,double *A,int m,int n){

  double *pdct;
  pdct=(double *) malloc(m*n*sizeof(double));
  int i;
  int j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      if(d>0){
        *(pdct+(i*n+j))=pow(d,*(A+(i*n+j)));
      }
      else{
        *(pdct+(i*n+j))=-1 * pow(abs(d),*(A+(i*n+j)));
      }
    }
  }
  return pdct;
}

//reference: http://cboard.cprogramming.com/cplusplus-programming/30001-determinant-calculation.html
double matrx_det(double *A,int n){
  int i, j, k;
  double **m;
  double det = 1;
  m = (double **) malloc(n*sizeof(double *));
  for ( i = 0; i < n; i++ ){
    m[i] = (double *) malloc(n*sizeof(double));
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ ){
      m[i][j] = *(A+(i*n+j));
    }
  }
  for ( k = 0; k < n; k++ ) {
    if ( m[k][k] == 0 ) {
      int ok = 0;
      for ( j = k; j < n; j++ ) {
        if (m[j][k] != 0 ){
          ok = 1;
	}
      }
      if (ok==0){
        return 0;
      }
      for ( i = k; i < n; i++ ){
        double tmp= m[i][j];
       	m[i][j]=m[i][k];
	m[i][k]=tmp;
      }
      det = -det;
    }
    det *= m[k][k];
    if ( k + 1 < n ) {
      for ( i = k + 1; i < n; i++ ) {
        for ( j = k + 1; j < n; j++ ){
          m[i][j] = m[i][j] - m[i][k] * m[k][j] / m[k][k];
        }
      }
    }
  }
  for ( i = 0; i < n; i++ ){
    free(m[i]);
  }
  free(m);
  return det;
}

double* ones(int m,int n){

  double *o;
  o = (double *) malloc(m*n*sizeof(double));
  int i;
  for(i=0;i<m*n;i++){
    *(o+i)=1;
  }
  return o;
}

double* eye(int n){

  double *iden;
  iden = (double *) malloc(n * n * sizeof(double));
  memset(iden,0,n*n*sizeof(double));
  int i;
  for(i=0;i<n;i++){
        *(iden+(i*n+i))=1;
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

double* tran(double *A,int m,int n){

  double* t;
  t = (double *) malloc(m*n*sizeof(double));
  int i;
  int j;
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      *(t+(i+j*m))=*(A+(i*n+j));
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
