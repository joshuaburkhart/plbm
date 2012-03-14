#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "jlapack.h"
#include "asa047.h"

extern double a_initVh[];
extern double a_initVp[];
extern double a_X[];
extern double a_tau1[];
extern double a_tau2[];
extern double n;
extern double p;
extern double q;

double funct(double *d1_d2);

int main(int argc,char *argv[]) {
    double *d1_d2;
    int i;
    d1_d2 = (double *) malloc(2 * sizeof(double));
    if(argc!=3) {
        printf("argument mismatch: you must supply d1 and d2 as doubles...aborting\n\n");
        return -1;
    }
    if(sizeof(double)==sizeof(argv[1])) {
        *(d1_d2) = atof(argv[1]);
    } else {
        return -2;
    }
    if(sizeof(double)==sizeof(argv[2])) {
        *(d1_d2+1) = atof(argv[2]);
    } else {
        return -3;
    }

    printf("funct returns %f\n",funct(d1_d2));
    
return 0;

    //reference: http://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html
    //reference: www.scholarpedia.org/article/Nelder-Mead_algorithm
    //reference: http://tolstoy.newcastle.edu.au/R/help/06/06/28963.html
    double XMIN[2]; //coordinates of minimum value
    double YNEWLO; //minimum value
    double REQMIN = 0.0001; //termination variance limit
    int KONVGE = 100; //frequency of convergence tests
    int KCOUNT = 10000; //max number of iterations
    int ICOUNT; //number of evaluations
    int NUMRES; //number of restarts
    int IFAULT; //error indicator
    nelmin(funct,2,d1_d2,XMIN,&YNEWLO,REQMIN,d1_d2,KONVGE,KCOUNT,&ICOUNT,&NUMRES,&IFAULT);
   
    printf("finished nelmin\n");

    printf("minimization coordinates: %f %f\n",*(XMIN),*(XMIN+1));
    printf("minimum value: %f\n",YNEWLO);

    free(d1_d2);
    return 0;
}

double funct(double *d1_d2) {
    
    double d1 = *(d1_d2);
    double d2 = *(d1_d2+1);

    /*
    output(a_initVh,p,p);
    output(a_initVp,q,q);
    printf("n: %f\n",n);
    printf("p: %f\n",p);
    printf("q: %f\n",q);
    output(a_X,n,1);
    output(a_tau1,p,p);
    output(a_tau2,q,q);
    printf("d1: %f\nd2: %f\n",d1,d2);
    */

//Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);--------------------------------------Vh

    double *A = array_pow(d1,a_tau1,p,p);
    double *B = matrx_mlt(2,a_initVh,p,p);
    double *C = array_pow(d1,B,p,p);
    double *D = matrx_sub(1,C,p,p);
    double *E = array_mlt(A,p,p,D);
    double d1sq = d1*d1;
    double omd1 = 1 - d1sq;
    double *Vh = array_rdv(E,p,p,omd1);

    output(A,p,p);
    output(B,p,p);
    output(C,p,p);
    output(D,p,p);
    output(E,p,p);
    printf("d1sq: %f\n",d1sq);
    printf("omd1: %f\n",omd1);
    output(Vh,p,p);

//Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);-------------------------------------Vp

    A = array_pow(d2,a_tau2,q,q);
    B = matrx_mlt(2,a_initVp,q,q);
    C = array_pow(d2,B,q,q);
    D = matrx_sub(1,C,q,q);
    E = array_mlt(A,q,q,D);
    double d2sq = d2*d2;
    double omd2 = 1 - d2sq;
    double *Vp = array_rdv(E,q,q,omd2);

    output(A,q,q);
    output(B,q,q);
    output(C,q,q);
    output(D,q,q);
    output(E,q,q);
    printf("d2sq: %f\n",d2sq);
    printf("omd2: %f\n",omd2);
    output(Vp,q,q);

//Vh=Vh./det(Vh)^(1/p);-----------------------------------Vh
    double dtm = matrx_det(Vh,p);
    printf("dtm = %f\n",dtm);
    double dtm21op = pow(dtm,(1/p));
    Vh = array_rdv(Vh,p,p,dtm21op);

    /*
    output(Vh,p,p);
    */

//Vp=Vp./det(Vp)^(1/q);-----------------------------------Vp

    dtm = matrx_det(Vp,q);
    printf("dtm = %f\n",dtm);
    double dtm21oq = pow(dtm,(1/q));
    Vp = array_rdv(Vp,q,q,dtm21oq);

    /*
    output(Vp,q,q);
    */

//V=kron(Vp,Vh);-----------------------------------V

    double *V = kron(Vp,q,q,Vh,p,p);

    /*
    output(V,n,n);
    */

//invV=V\eye(n);-----------------------------------invV
//reference: http://www.netlib.org/lapack/double

    A = tran(V,n,n); //row major -> column major
    double *invV;
    int N=n;
    int lda=N;
    int ipiv[N];
    int info;
    int lwork=N*N;
    double work[lwork];
    dgetrf_(&N,&N,A,&lda,ipiv,&info);
    if(info!=0) {
        printf("dgetrf returns info code %i\n",info);
    }
    dgetri_(&N,A,&lda,ipiv,work,&lwork,&info);
    if(info!=0) {
        printf("dgetri returns info code %i\n",info);
    }

    /*
    int nrhs=N;
    B = eye(n);
    int ldb=N;
    dgesv_(&N,&nrhs,A,&lda,ipiv,B,&ldb,&info);
    printf("dgesv returns info status %i\n",info);
    */

    invV = tran(A,n,n); //column major -> row major

    /*
    output(invV,n,n);
    */

//U=ones(length(X),1);-----------------------------------U

    double *U = ones(n,1);

    /*
    output(U,n,1);
    */

//b=(U'*invV*U)\(U'*invV*X);-----------------------------------b

    A = tran(U,n,1);
    B = matrx_mlt2(A,1,n,invV,n,n);
    C = matrx_mlt2(B,1,n,U,n,1);
    D = matrx_mlt2(B,1,n,a_X,n,1);

    double c = *(C);//C should be a 1 x 1 matrix
    double d = *(D);//D should be a 1 x 1 matrix
    double b = d/c;

    
    printf("c = %f\n",c);
    printf("d = %f\n",d);
    printf("b = %f\n",b);
    

//H=X-b;-----------------------------------H

    double *H = matrx_sub3(a_X,n,1,b);

    /*
    output(H,n,1);
    */

//MSE=(H'*invV*H)/(n-1);-----------------------------------MSE

    A = tran(H,n,1);
    B = matrx_mlt2(A,1,n,invV,n,n);
    C = matrx_mlt2(B,1,n,H,n,1);
    c = *(C);//C should be a 1 x 1 matrix
    double MSE = c/(n -1);

    /*
    printf("MSE = %f\n",MSE);
    */

    //freeing up space
    free(H);
    free(invV);
    free(Vp);
    free(V);
    free(A);
    free(B);
    free(C);
    free(D);
    free(E);
    free(Vh);
    return MSE;
}
