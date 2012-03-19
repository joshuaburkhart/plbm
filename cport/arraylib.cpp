#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "./lib/lapack.h"
#include "./lib/arraylib.h"

void matrx_inv(double out[],double *A,int n) {

    tran(out,A,n,n);

    ptrdiff_t N=n;
    ptrdiff_t M=n;
    ptrdiff_t lda=n;
    ptrdiff_t ipiv[N];
    ptrdiff_t info=0;
    ptrdiff_t lwork=n*n;
    double work[lwork];

    dgetrf(&M,&N,out,&lda,ipiv,&info);
    if(info!=0) {
        //printf("dgetrf returns info code %i ... inverse could not be calculated\n",info);
        //printf("M:   %i\n",M);
        //printf("N:   %i\n",N);
        //printf("lda: %i\n",lda);
        info=0;
    }

    dgetri(&N,out,&lda,ipiv,work,&lwork,&info);
    if(info!=0) {
        //printf("dgetri returns info code %i ... inverse could not be calculated\n",info);
        //printf("N:   %i\n",N);
        //printf("lda: %i\n",lda);
    }

    tran(out,out,n,n);
    return;
}

void array_rdv(double out[],double *A,int m,int n,double d) {

    #pragma omp parallel for
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            out[i*n+j] =*(A+(i*n+j)) / d;
        }
    }
    return;
}

void matrx_sub3(double out[],double *A,int m,int n,double d) {

    #pragma omp parallel for
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            out[i*n+j] =*(A+(i*n+j)) - d;
        }
    }
    return;
}

void matrx_sub(double out[],double d,double *A,int m,int n) {

    #pragma omp parallel for
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            out[i*n+j] =d - *(A+(i*n+j));
        }
    }
    return;
}

void matrx_mlt2(double out[],double *A,int ma,int na,double *B,int mb,int nb) {

    #pragma omp parallel for
    for(int i=0; i<ma; i++) {
        for(int j=0; j<nb; j++) {
            double sum=0;
            for(int k=0; k<na; k++) {
                sum+=*(A+(i*na+k)) * *(B+(j+k*nb));
            }
            out[i*nb+j] =sum;
        }
    }
    return;
}

void matrx_mlt(double out[],double d,double *A,int m,int n) {

    #pragma omp parallel for
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            out[i*n+j]  = *(A+(i*n+j)) * d;
        }
    }
    return;
}

void array_mlt(double out[],double A[],int m,int n,double *B) {

    #pragma omp parallel for
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            out[i*n+j] =A[i*n+j] * *(B+(i*n+j));
        }
    }
    return;
}

void array_pow(double out[],double d,double *A,int m,int n) {

    #pragma omp parallel for
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            out[i*n+j] =pow(d,*(A+(i*n+j)));
        }
    }
    return;
}

double matrx_det(double *A,int n) {
    double luout[n*n];
    tran(luout,A,n,n);
    ptrdiff_t N=n;
    ptrdiff_t M=n;
    ptrdiff_t lda=n;
    ptrdiff_t ipiv[N];
    ptrdiff_t info=0;
    ptrdiff_t lwork=N*N;
    double work[lwork];
    dgetrf(&M,&N,luout,&lda,ipiv,&info);
    if(info!=0) {
        //printf("dgetrf returns info code %i ... determinant cannot be calculated\n",info);
        //printf("M:   %i\n",M);
        //printf("N:   %i\n",N);
        //printf("lda: %i\n",lda);
    }
    double luout2[n*n];
    tran(luout2,luout,n,n);
    double diag=1;
    for(int i = 0; i < n; i++) {
        diag *= luout2[i * n + i];
    }
    #pragma omp parallel for
    for(int i = 0; i < n; i++) {
        luout2[i *  n + i]  = 1;
        for(int j = i+1; j < n; j++) {
            luout2[i *  n + j]  = 0;
        }
    }
    double dtm=det_l(luout2,n);
    return(dtm * diag);
}

static double det_l(double *A,int n) {
    int i, j, k;
    double m[n][n];
    double det = 1;
    #pragma omp parallel for
    for (int x = 0; x < n; x++ ) {
        for (int y = 0; y < n; y++ ) {
            m[x][y] = *(A+(x*n+y));
        }
    }
    for ( k = 0; k < n; k++ ) {
        if ( m[k][k] == 0 ) {
            int ok = 0;
            for ( j = k; j < n; j++ ) {
                if (m[j][k] != 0 ) {
                    ok = 1;
                }
            }
            if (ok==0) {
                return 0;
            }
            for ( i = k; i < n; i++ ) {
                double tmp= m[i][j];
                m[i][j]=m[i][k];
                m[i][k]=tmp;
            }
            det = -det;
        }
        det *= m[k][k];
        if ( k + 1 < n ) {
            for ( i = k + 1; i < n; i++ ) {
                for ( j = k + 1; j < n; j++ ) {
                    m[i][j] = m[i][j] - m[i][k] * m[k][j] / m[k][k];
                }
            }
        }
    }
    return det;
}

void ones(double out[],int m,int n) {

    #pragma omp parallel for
    for(int i=0; i<m*n; i++) {
        out[i]=1;
    }
    return;
}

void kron(double out[],double *A,int ma,int na,double *B,int mb,int nb) {

    int min_m = ma < mb ? ma : mb;
    int min_n = na < nb ? na : nb;
    #pragma omp parallel for
    for(int i=0; i<ma*mb; i++) {
        for(int j=0; j<na*nb; j++) {
            int a_row = i/min_m;
            int a_col = j/min_n;
            int b_row = i%min_m;
            int b_col = j%min_n;
            double val=*(A+(a_row*na+a_col)) * *(B+(b_row*nb+b_col));
            out[i*na*nb +j]=val;
        }
    }
    return;
}

void tran(double out[],double *A,int m,int n) {

    #pragma omp parallel for
    for(int i=0; i<m; i++) {
        for(int j=0; j<n; j++) {
            out[i+j*m] =*(A+(i*n+j));
        }
    }
    return;
}
