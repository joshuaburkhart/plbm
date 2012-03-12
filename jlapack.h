#define L 4

extern double a_initVh[];
extern double a_initVp[];
extern double n;
extern double p;
extern double q;
extern double a_X[];
extern double a_tau1[];
extern double a_tau2[];

void output(double *matrix,int m,int n);
double* eye(int n);
double* kron(double *A,int ma,int na,double *B,int mb,int nb);
double* tran(double *A,int m,int n);
double* ones(int m,int n);
double matrx_det(double *A,int n);
//reference: http://www.netlib.org/lapack/double/dgesv.f
void dgesv_(const int *N, const int *nrhs,double *A,const int *lda,int *ipiv,double *b,const int *ldb,int *info);
double* array_pow(double d,double *A,int m,int n);
double* array_mlt(double *A,int m,int n,double *B);
double* matrx_mlt(double d,double *A,int ma,int na);
double* matrx_mlt2(double *A,int ma,int na,double *B,int mb,int nb);
double* matrx_sub(double d,double *A,int m,int n);
double* matrx_sub2(double *A,int m,int n,double *B);
double* array_rdv(double *A,int m,int n,double d);
