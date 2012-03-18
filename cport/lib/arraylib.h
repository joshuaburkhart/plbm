
double matrx_det(double *A,int n);
static double det_l(double *A,int n);

//double* matrx_inv(double *A,int n);
void matrx_inv(double out[],double *A,int n);
//double* eye(int n);
//double* kron(double *A,int ma,int na,double *B,int mb,int nb);
void kron(double out[],double *A,int ma,int na,double *B,int mb,int nb);
//double* tran(double *A,int m,int n);
void tran(double out[],double *A,int m,int n);
//double* ones(int m,int n);
void ones(double out[],int m,int n);
//double* array_pow(double d,double *A,int m,int n);
void array_pow(double out[],double d,double *A,int m,int n);
//double* array_mlt(double *A,int m,int n,double *B);
void array_mlt(double out[],double A[],int m,int n,double *B);
//double* matrx_mlt(double d,double *A,int ma,int na);
void matrx_mlt(double out[],double d,double *A,int m,int n);
//double* matrx_mlt2(double *A,int ma,int na,double *B,int mb,int nb);
void matrx_mlt2(double out[],double *A,int ma,int na,double *B,int mb,int nb);
//double* matrx_sub(double d,double *A,int m,int n);
void matrx_sub(double out[],double d,double *A,int m,int n);
//double* matrx_sub2(double *A,int m,int n,double *B);
//double* matrx_sub3(double *A,int m,int n,double d);
void matrx_sub3(double out[],double *A,int m,int n,double d);
void array_rdv(double out[],double *A,int m,int n,double d);

/*
void dgesv_(const int *N, const int *nrhs,double *A,const int *lda,int *ipiv,double *b,const int *ldb,int *info);
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetrf(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
*/
