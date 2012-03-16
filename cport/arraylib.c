#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "./lib/arraylib.h"

double* array_rdv(double *A,int m,int n,double d) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=*(A+(i*n+j)) / d;
        }
    }
    return diff;
}

double* matrx_sub3(double *A,int m,int n,double d) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=*(A+(i*n+j)) - d;
        }
    }
    return diff;
}


double* matrx_sub2(double *A,int m,int n,double *B) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=*(A+(i*n+j)) - *(B+(i*n+j));
        }
    }
    return diff;
}

double* matrx_sub(double d,double *A,int m,int n) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=d - *(A+(i*n+j));
        }
    }
    return diff;
}

double* matrx_mlt2(double *A,int ma,int na,double *B,int mb,int nb) {

    double *pdct;
    pdct=(double *) malloc(ma*nb*sizeof(double));
    int i;
    int j;
    int k;
    for(i=0; i<ma; i++) {
        for(j=0; j<nb; j++) {
            double sum=0;
            for(k=0; k<na; k++) {
                sum+=*(A+(i*na+k)) * *(B+(j+k*nb));
            }
            *(pdct+(i*nb+j))=sum;
        }
    }
    return pdct;
}

double* matrx_mlt(double d,double *A,int m,int n) {

    double *pdct;
    pdct=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(pdct+(i*n+j))=*(A+(i*n+j)) * d;
        }
    }
    return pdct;
}

double* array_mlt(double *A,int m,int n,double *B) {

    double *pdct;
    pdct=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(pdct+(i*n+j))=*(A+(i*n+j)) * *(B+(i*n+j));
        }
    }
    return pdct;
}

double* array_pow(double d,double *A,int m,int n) {

    double *pdct;
    pdct=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(pdct+(i*n+j))=pow(d,*(A+(i*n+j)));
        }
    }
    return pdct;
}

double matrx_det(double *A,int n) {
    double *lu;
    lu = malloc(n*n*sizeof(double));
    lu = tran(A,n,n); /*row major -> column major*/
    int N=n;
    int lda=n;
    int ipiv[N];
    int info;
    int lwork=N*N;
    double work[lwork];
    dgetrf_(&N,&N,lu,&lda,ipiv,&info);
    if(info!=0) {
        printf("dgetrf returns info code %i ... determinant cannot be calculated\n",info);
        printf("N:   %i\n",N);
	printf("lda: %i\n",lda);
    }
    lu = tran(lu,n,n);
    double diag=1;
    int i;
    for(i = 0; i < n; i++) {
        double multiplier =  *(lu + (i * n + i));
        diag *= multiplier;
    }
    for(i = 0; i < n; i++) {
        *(lu + (i *  n + i)) = 1;
        int j;
        for(j = i+1; j < n; j++) {
            *(lu + ( i *  n + j)) = 0;
        }
    }
    double dtm=det_l(lu,n);
    free(lu);
    return(dtm * diag);
}

static double det_l(double *A,int n) {
    int i, j, k;
    double **m;
    double det = 1;
    m = (double **) malloc(n*sizeof(double *));
    for ( i = 0; i < n; i++ ) {
        m[i] = (double *) malloc(n*sizeof(double));
    }
    for ( i = 0; i < n; i++ ) {
        for ( j = 0; j < n; j++ ) {
            m[i][j] = *(A+(i*n+j));
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
    for ( i = 0; i < n; i++ ) {
        free(m[i]);
    }
    free(m);
    return det;
}

double* ones(int m,int n) {

    double *o;
    o = (double *) malloc(m*n*sizeof(double));
    int i;
    for(i=0; i<m*n; i++) {
        *(o+i)=1;
    }
    return o;
}

double* eye(int n) {

    double *iden;
    iden = (double *) malloc(n * n * sizeof(double));
    memset(iden,0,n*n*sizeof(double));
    int i;
    for(i=0; i<n; i++) {
        *(iden+(i*n+i))=1;
    }
    return iden;
}

double* kron(double *A,int ma,int na,double *B,int mb,int nb) {

    double* k;
    k = (double *) malloc(ma*mb*na*nb*sizeof(double));
    int i;
    int j;
    int min_m = ma < mb ? ma : mb;
    int min_n = na < nb ? na : nb;
    for(i=0; i<ma*mb; i++) {
        for(j=0; j<na*nb; j++) {
            int a_row = i/min_m;
            int a_col = j/min_n;
            int b_row = i%min_m;
            int b_col = j%min_n;
            double val=*(A+(a_row*na+a_col)) * *(B+(b_row*nb+b_col));
            *(k+((i*na*nb)+j))=val;
        }
    }
    return k;
}

double* tran(double *A,int m,int n) {

    double* t;
    t = (double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(t+(i+j*m))=*(A+(i*n+j));
        }
    }
    return t;
}

void output(double *matrix,int m,int n) {

    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            printf("%f ",*(matrix+(i*n+j)));
        }
        printf("\n");
    }
}
