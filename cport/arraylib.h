
void output(double *matrix,int m,int n);
double matrx_det(double *A,int n);
static double det_l(double *A,int n);

double* eye(int n);
double* kron(double *A,int ma,int na,double *B,int mb,int nb);
double* tran(double *A,int m,int n);
double* ones(int m,int n);
double* array_pow(double d,double *A,int m,int n);
double* array_mlt(double *A,int m,int n,double *B);
double* matrx_mlt(double d,double *A,int ma,int na);
double* matrx_mlt2(double *A,int ma,int na,double *B,int mb,int nb);
double* matrx_sub(double d,double *A,int m,int n);
double* matrx_sub2(double *A,int m,int n,double *B);
double* matrx_sub3(double *A,int m,int n,double d);
double* array_rdv(double *A,int m,int n,double d);

void dgesv_(const int *N, const int *nrhs,double *A,const int *lda,int *ipiv,double *b,const int *ldb,int *info);
void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

