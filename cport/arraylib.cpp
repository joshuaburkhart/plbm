/*
   distributed under the terms of the GNU General Public License
   Copyright 2012 Joshua Burkhart
   */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "./lib/lapack.h"
#include "./lib/arraylib.h"

void matrx_inv(double out[],double *A,int n) {

	ptrdiff_t N=n;
	ptrdiff_t M=n;
	ptrdiff_t lda=n;
	ptrdiff_t ipiv[N];
	ptrdiff_t info=0;
	ptrdiff_t lwork=n*n;
	double work[lwork];

	dgetrf(&M,&N,out,&lda,ipiv,&info);
	dgetri(&N,out,&lda,ipiv,work,&lwork,&info);

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
				sum += *(A+(i*na+k)) * *(B+(j+k*nb));
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

void array_mlt(double out[],double A[],int m,int n,double B[]) {

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

	double B[n*n];

#pragma omp parallel for
	for(int i=0; i < n*n; i++){
		B[i] = *(A+i);
	}

	double mult;
	double deter=1.0;
	for(int i=0;i<n;i++){   
		for(int j=0;j<n;j++){   
			mult= *(B+(j*n +i)) / *(B+(i*n +i));
			if(i!=j){
#pragma omp parallel for
				for(int k=0;k<n;k++){   
					*(B+(j*n +k)) = *(B+(j*n +k)) - (*(B+(i*n +k)) * mult);
				}
			}	
		}
	}

	for(int i=0;i<n;i++){
		deter *= *(B+(i*n +i));
	}   
	return(deter);
}

void ones(double out[],int m,int n) {

#pragma omp parallel for
	for(int i=0; i<m*n; i++) {
		out[i]=1.0;
	}
	return;
}

void kron(double out[],double *A,int ma,int na,double *B,int mb,int nb) {

	int min_m = ma < mb ? ma : mb;
	int min_n = na < nb ? na : nb;

#pragma omp parallel for
	for(int i=0; i<ma*mb; i++) {
		for(int j=0; j<na*nb; j++) {
			out[i*na*nb +j]=*(A+(i/min_m*na+j/min_n)) * *(B+(i%min_m*nb+j%min_n));
		}
	}
	return;
}

void tran(double out[],double *A,int m,int n) {

	double tmp[m*n];

#pragma omp parallel for
	for(int i=0; i<m; i++) {
		for(int j=0; j<n; j++) {
			tmp[j*m + i] = *(A+(i*n + j));
		}
	}

#pragma omp parallel for
	for(int i=0;i <(m*n); i++){
		out[i] = tmp[i];
	}
	return;
}
