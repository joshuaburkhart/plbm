/*
distributed under the terms of the GNU General Public License
Copyright 2012 Joshua Burkhart
*/

double matrx_det(double *A,int n);
static double det_l(double *A,int n);

void matrx_inv(double out[],double *A,int n);
void kron(double out[],double *A,int ma,int na,double *B,int mb,int nb);
void tran(double out[],double *A,int m,int n);
void ones(double out[],int m,int n);
void array_pow(double out[],double d,double *A,int m,int n);
void array_mlt(double out[],double A[],int m,int n,double B[]);
void matrx_mlt(double out[],double d,double *A,int m,int n);
void matrx_mlt2(double out[],double *A,int ma,int na,double *B,int mb,int nb);
void matrx_sub(double out[],double d,double *A,int m,int n);
void matrx_sub3(double out[],double *A,int m,int n,double d);
void array_rdv(double out[],double *A,int m,int n,double d);
