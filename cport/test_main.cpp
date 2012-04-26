/*
distributed under the terms of the GNU General Public License
Copyright 2012 Joshua Burkhart
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./lib/lapack.h"
#include "./lib/arraylib.h"
#include "./lib/minimize.h"

extern double *initVh;
extern double *initVp;
extern double *X;
extern double *tau1;
extern double *tau2;
extern int n;
extern int p;
extern int q;

void output(double *matrix,int m,int n);
void test_funct(double d1, double d2, double correct_return, double epsilon);

int main(int argc,char *argv[]) {
    //test each arraylib function with several values

    /*Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);--------------------------------------*/

    double A_PP[p*p];
    matrx_mlt(A_PP,2.00,initVh,p,p);
    array_pow(A_PP,d1,A_PP,p,p);
    matrx_sub(A_PP,1.00,A_PP,p,p);

    double B_PP[p*p];
    array_pow(B_PP,d1,tau1,p,p);

    array_mlt(B_PP,B_PP,p,p,A_PP);
    array_rdv(B_PP,B_PP,p,p,1.00 - d1*d1);

    double output = {1.1776,0.7765,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0.7765,1.1776,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,  0.4960  ,  0.0670
                     0 ,        0  ,  1.1776 ,        0  ,       0 ,        0 ,        0  ,       0 ,        0  ,       0  ,       0 ,        0
                     0.1629  ,  0.1629  ,       0 ,   1.1776  ,  0.1629 ,   0.3000 ,   0.1629  ,  0.0670 ,   0.0670  ,  0.0670  ,  0.1629 ,   0.0670
                     0.3000   , 0.3000  ,       0 ,   0.1629  ,  1.1776 ,   0.1629 ,   0.3000  ,  0.0670 ,   0.0670  ,  0.0670  ,  0.3000 ,   0.0670
                     0.1629   , 0.1629  ,       0 ,   0.3000  ,  0.1629 ,   1.1776 ,   0.1629  ,  0.0670 ,   0.0670  ,  0.0670  ,  0.1629 ,   0.0670
                     0.4960   , 0.4960  ,       0 ,   0.1629  ,  0.3000 ,   0.1629 ,   1.1776  ,  0.0670 ,   0.0670  ,  0.0670  ,  0.7765 ,   0.0670
                     0.0670   , 0.0670  ,       0 ,   0.0670  ,  0.0670 ,   0.0670 ,   0.0670  ,  1.1776 ,   0.3000  ,  0.3000  ,  0.0670 ,   0.1629
                     0.0670   , 0.0670  ,       0 ,   0.0670  ,  0.0670 ,   0.0670 ,   0.0670  ,  0.3000 ,   1.1776  ,  0.4960  ,  0.0670 ,   0.1629
                     0.0670   , 0.0670  ,       0 ,   0.0670  ,  0.0670 ,   0.0670 ,   0.0670  ,  0.3000 ,   0.4960  ,  1.1776  ,  0.0670 ,   0.1629
                     0.4960   , 0.4960  ,       0 ,   0.1629  ,  0.3000 ,   0.1629 ,   0.7765  ,  0.0670 ,   0.0670  ,  0.0670  ,  1.1776 ,   0.0670
                     0.0670   , 0.0670  ,       0 ,   0.0670  ,  0.0670 ,   0.0670 ,   0.0670  ,  0.1629 ,   0.1629  ,  0.1629  ,  0.0670 ,   1.1776
                    }

                    /*Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);-------------------------------------*/
                    /*
                        double A_QQ[q*q];
                        matrx_mlt(A_QQ,2.00,initVp,q,q);
                        array_pow(A_QQ,d2,A_QQ,q,q);
                        matrx_sub(A_QQ,1.00,A_QQ,q,q);

                        double B_QQ[q*q];
                        array_pow(B_QQ,d2,tau2,q,q);

                        array_mlt(B_QQ,B_QQ,q,q,A_QQ);
                        array_rdv(B_QQ,B_QQ,q,q,1.00 - d2*d2);
                    */
                    /*Vh=Vh./det(Vh)^(1/p);-----------------------------------*/
                    /*
                        array_rdv(B_PP,B_PP,p,p,pow(matrx_det(B_PP,p),(1.00/((double) p))));
                    */
                    /*Vp=Vp./det(Vp)^(1/q);-----------------------------------*/
                    /*
                        array_rdv(B_QQ,B_QQ,q,q,pow(matrx_det(B_QQ,q),(1.00/((double) q))));
                    */
                    /*V=kron(Vp,Vh);-----------------------------------*/
                    /*
                        double *V_NN;
                        V_NN = (double *) malloc(n*n*sizeof(double));
                        kron(V_NN,B_QQ,q,q,B_PP,p,p);
                    */
                    /*invV=V\eye(n);-----------------------------------*/
                    /*
                        matrx_inv(V_NN,V_NN,n); //saving memory by reusing name
                    */
                    /*U=ones(length(X),1);-----------------------------------*/
                    /*
                        double A_N[n];
                        ones(A_N,n,1.00);
                    */
                    /*b=(U'*invV*U)\(U'*invV*X);-----------------------------------*/
                    /*
                        double B_N[n];
                        double B_NN[n*n];
                        tran(B_N,A_N,n,1.00);
                        matrx_mlt2(B_NN,B_N,1.00,n,V_NN,n,n);
                        matrx_mlt2(B_N,B_NN,1.00,n,X,n,1.00);
                        matrx_mlt2(A_N,B_NN,1.00,n,A_N,n,1.00);
                    */
                    /*H=X-b;-----------------------------------*/
                    /*
                        matrx_sub3(B_N,X,n,1.00,B_N[0] / A_N[0]);
                    */
                    /*MSE=(H'*invV*H)/(n-1);-----------------------------------*/
                    /*
                        tran(A_N,B_N,n,1.00);
                        matrx_mlt2(B_NN,A_N,1.00,n,V_NN,n,n);
                        matrx_mlt2(B_NN,B_NN,1.00,n,B_N,n,1.00);

                        free(V_NN);
                        return B_NN[0] / ((double) n - 1);
                    */

                    //test funct with several values

                    double epsilon=0.00000001;

    test_funct(0.5, 0.5, 0.0101290635, epsilon);
    test_funct(0.0, 0.5, 0.0116081347, epsilon);
    test_funct(0.5, 0.0, 0.0098632064, epsilon);
    test_funct(1.0, 0.5, NAN, epsilon);
    test_funct(0.5, 1.0, NAN, epsilon);
    test_funct(0.0, 0.0, 0.0109859425, epsilon);
    test_funct(1.0, 1.0, NAN, epsilon);
    test_funct(0.75, 0.75, 0.0106599215, epsilon);
    test_funct(0.25, 0.25, 0.0100825613, epsilon);

    //test nelmin with several values

    double *d1_d2;
    d1_d2 = (double *) malloc(2 * sizeof(double));

    *(d1_d2) = .5;
    *(d1_d2+1) = .5;

    /*TODO: play with konvge settings for optimization*/
    double STEP[2];
    STEP[0] = (*(d1_d2 )  == 0 ? 0.00025 : 0.95 * *(d1_d2));
    STEP[1] = (*(d1_d2+1) == 0 ? 0.00025 : 0.95 * *(d1_d2+1));
    double XMIN[2];         /*coordinates of minimum value*/
    double YNEWLO;          /*minimum value*/
    double REQMIN = 0.0001; /*termination variance limit*/
    int KONVGE = 100;       /*frequency of convergence tests*/
    int KCOUNT = 10000;     /*max number of iterations*/
    int ICOUNT;             /*number of evaluations*/
    int NUMRES;             /*number of restarts*/
    int IFAULT;             /*error indicator*/

    nelmin(funct,2,d1_d2,XMIN,&YNEWLO,REQMIN,STEP,KONVGE,KCOUNT,&ICOUNT,&NUMRES,&IFAULT);

    /* /////////////////////// */
    /* print results to screen */
    /* /////////////////////// */

    printf("est=\n  %f  %f\n",*(XMIN),*(XMIN+1));
    printf("MSE=\n  %f\n",YNEWLO);

    free(d1_d2);
    return 0;
}

void test_funct(double d1, double d2, double correct_return, double epsilon) {

    double *d1_d2;
    d1_d2 = (double *) malloc(2 * sizeof(double));

    *(d1_d2) = d1;
    *(d1_d2+1) = d2;

    double mse=funct(d1_d2);
    double error = mse - correct_return;
    if( isnan(mse) && isnan(correct_return) || (error * error) < epsilon ) {
        printf("* PASSED");
    } else {
        printf("X FAILED");
    }
    printf(" -> funct()\n");
    printf("\t%f returned\n",mse);
    printf("\t%f correct\n",correct_return);

    free(d1_d2);
}

void test_vh_calculation(double d1, double d2, double *correct_return[], double epsilon) {

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
