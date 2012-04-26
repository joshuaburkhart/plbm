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

int all_passed = 1;

void output(double *matrix,int m,int n);
void test_funct(double d1, double d2, double correct_return, double epsilon);
void test_vh_calculation(double d1, double d2, double correct_return[], int correct_size, double epsilon);

int main(int argc,char *argv[]) {

    printf("\n\ttesting...\n\n");
    double epsilon=0.00000001;

    //test each arraylib function with several values

    /*Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);--------------------------------------*/

    double expect_vh_5_5[]= {1.1776,0.7765,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0.7765,1.1776,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0,0,1.1776,0,0,0,0,0,0,0,0,0,0.1629,0.1629,0,1.1776,0.1629 ,0.3000 ,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.3000,0.3000,0,0.1629,1.1776 ,0.1629 ,0.3000,0.0670,0.0670,0.0670,0.3000,0.0670,0.1629,0.1629,0,0.3000,0.1629 ,1.1776 ,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.4960,0.4960,0,0.1629,0.3000 ,0.1629 ,1.1776,0.0670,0.0670,0.0670,0.7765,0.0670,0.0670,0.0670,0,0.0670,0.0670 ,0.0670 ,0.0670,1.1776,0.3000,0.3000,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670 ,0.0670 ,0.0670,0.3000,1.1776,0.4960,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670 ,0.0670 ,0.0670,0.3000,0.4960,1.1776,0.0670,0.1629,0.4960,0.4960,0,0.1629,0.3000 ,0.1629 ,0.7765,0.0670,0.0670,0.0670,1.1776,0.0670,0.0670,0.0670,0,0.0670,0.0670 ,0.0670 ,0.0670,0.1629,0.1629,0.1629,0.0670,1.1776};
    int expect_vh_5_5_size = sizeof(expect_vh_5_5) / sizeof(double);
    test_vh_calculation(0.5,0.5,expect_vh_5_5,expect_vh_5_5_size,epsilon);
    
    double expect_vh_0_5[]= {1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1};
    int expect_vh_0_5_size = sizeof(expect_vh_0_5) / sizeof(double);
    test_vh_calculation(0.0,0.5,expect_vh_0_5,expect_vh_0_5_size,epsilon);
    
    double expect_vh_5_0[]= {1.1776,0.7765,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0.7765,1.1776,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0,0,1.1776,0,0,0,0,0,0,0,0,0,0.1629,0.1629,0,1.1776,0.1629,0.3000,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.3000,0.3000,0,0.1629,1.1776,0.1629,0.3000,0.0670,0.0670,0.0670,0.3000,0.0670,0.1629,0.1629,0,0.3000,0.1629,1.1776,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.4960,0.4960,0,0.1629,0.3000,0.1629,1.1776,0.0670,0.0670,0.0670,0.7765,0.0670,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,1.1776,0.3000,0.3000,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.3000,1.1776,0.4960,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.3000,0.4960,1.1776,0.0670,0.1629,0.4960,0.4960,0,0.1629,0.3000,0.1629,0.7765,0.0670,0.0670,0.0670,1.1776,0.0670,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.1629,0.1629,0.1629,0.0670,1.1776};
    int expect_vh_5_0_size = sizeof(expect_vh_5_0) / sizeof(double);
    test_vh_calculation(0.5,0.0,expect_vh_5_0,expect_vh_5_0_size,epsilon);
    
    double expect_vh_1_5[]= {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
    int expect_vh_1_5_size = sizeof(expect_vh_1_5) / sizeof(double);
    test_vh_calculation(1.0,0.5,expect_vh_1_5,expect_vh_1_5_size,epsilon);

    double expect_vh_5_1[]= {1.1776,0.7765,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0.7765,1.1776,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0,0,1.1776,0,0,0,0,0,0,0,0,0,0.1629,0.1629,0,1.1776,0.1629,0.3000,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.3000,0.3000,0,0.1629,1.1776,0.1629,0.3000,0.0670,0.0670,0.0670,0.3000,0.0670,0.1629,0.1629,0,0.3000,0.1629,1.1776,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.4960,0.4960,0,0.1629,0.3000,0.1629,1.1776,0.0670,0.0670,0.0670,0.7765,0.0670,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,1.1776,0.3000,0.3000,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.3000,1.1776,0.4960,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.3000,0.4960,1.1776,0.0670,0.1629,0.4960,0.4960,0,0.1629,0.3000,0.1629,0.7765,0.0670,0.0670,0.0670,1.1776,0.0670,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.1629,0.1629,0.1629,0.0670,1.1776};
    int expect_vh_5_1_size = sizeof(expect_vh_5_1) / sizeof(double);
    test_vh_calculation(0.5,1.0,expect_vh_5_1,expect_vh_5_1_size,epsilon);

    double expect_vh_0_0[]= {1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1};
    int expect_vh_0_0_size = sizeof(expect_vh_0_0) / sizeof(double);
    test_vh_calculation(0.0,0.0,expect_vh_0_0,expect_vh_0_0_size,epsilon);

    double expect_vh_1_1[]= {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
    int expect_vh_1_1_size = sizeof(expect_vh_1_1) / sizeof(double);
    test_vh_calculation(1.0,1.0,expect_vh_1_1,expect_vh_1_1_size,epsilon);

    double expect_vh_75_75[]= {1.3482,1.0327,0,0.3243,0.5263,0.3243,0.7607,0.1501,0.1501,0.1501,0.7607,0.1501,1.0327,1.3482,0,0.3243,0.5263,0.3243,0.7607,0.1501,0.1501,0.1501,0.7607,0.1501,0,0,1.3482,0,0,0,0,0,0,0,0,0,0.3243,0.3243,0,1.3482,0.3243,0.5263,0.3243,0.1501,0.1501,0.1501,0.3243,0.1501,0.5263,0.5263,0,0.3243,1.3482,0.3243,0.5263,0.1501,0.1501,0.1501,0.5263,0.1501,0.3243,0.3243,0,0.5263,0.3243,1.3482,0.3243,0.1501,0.1501,0.1501,0.3243,0.1501,0.7607,0.7607,0,0.3243,0.5263,0.3243,1.3482,0.1501,0.1501,0.1501,1.0327,0.1501,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,1.3482,0.5263,0.5263,0.1501,0.3243,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,0.5263,1.3482,0.7607,0.1501,0.3243,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,0.5263,0.7607,1.3482,0.1501,0.3243,0.7607,0.7607,0,0.3243,0.5263,0.3243,1.0327,0.1501,0.1501,0.1501,1.3482,0.1501,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,0.3243,0.3243,0.3243,0.1501,1.3482};
    int expect_vh_75_75_size = sizeof(expect_vh_75_75) / sizeof(double);
    test_vh_calculation(0.75,0.75,expect_vh_75_75,expect_vh_75_75_size,epsilon);

    double expect_vh_25_25[]= {1.0521,0.5069,0,0.0464,0.1100,0.0464,0.2403,0.0152,0.0152,0.0152,0.2403,0.0152,0.5069,1.0521,0,0.0464,0.1100,0.0464,0.2403,0.0152,0.0152,0.0152,0.2403,0.0152,0,0,1.0521,0,0,0,0,0,0,0,0,0,0.0464,0.0464,0,1.0521,0.0464,0.1100,0.0464,0.0152,0.0152,0.0152,0.0464,0.0152,0.1100,0.1100,0,0.0464,1.0521,0.0464,0.1100,0.0152,0.0152,0.0152,0.1100,0.0152,0.0464,0.0464,0,0.1100,0.0464,1.0521,0.0464,0.0152,0.0152,0.0152,0.0464,0.0152,0.2403,0.2403,0,0.0464,0.1100,0.0464,1.0521,0.0152,0.0152,0.0152,0.5069,0.0152,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,1.0521,0.1100,0.1100,0.0152,0.0464,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,0.1100,1.0521,0.2403,0.0152,0.0464,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,0.1100,0.2403,1.0521,0.0152,0.0464,0.2403,0.2403,0,0.0464,0.1100,0.0464,0.5069,0.0152,0.0152,0.0152,1.0521,0.0152,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,0.0464,0.0464,0.0464,0.0152,1.0521};
    int expect_vh_25_25_size = sizeof(expect_vh_25_25) / sizeof(double);
    test_vh_calculation(0.25,0.25,expect_vh_25_25,expect_vh_25_25_size,epsilon);

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
    if(all_passed == 1){
      printf("*** Congratulations! No tests FAILED, all PASSED!\n");
    } else {
      printf("XXX One or more tests FAILED\n");
    }

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
        all_passed = 0;
        printf("X FAILED");
    }
    printf(" -> funct()");
    printf(" -> d1 = %f",d1);
    printf(" -> d2 = %f",d2);
    printf(" -> %f returned",mse);
    printf(" -> %f correct\n",correct_return);

    free(d1_d2);
}

void test_vh_calculation(double d1, double d2, double correct_return[], int correct_size, double epsilon) {

    double A_PP[p*p];
    matrx_mlt(A_PP,2.00,initVh,p,p);
    array_pow(A_PP,d1,A_PP,p,p);
    matrx_sub(A_PP,1.00,A_PP,p,p);

    double B_PP[p*p];
    array_pow(B_PP,d1,tau1,p,p);

    array_mlt(B_PP,B_PP,p,p,A_PP);
    array_rdv(B_PP,B_PP,p,p,1.00 - d1*d1);
   
    double vh_size = sizeof(B_PP) / sizeof(double);

    if(vh_size == correct_size){
      for(int i=0; i <correct_size; i++){
        double actual = B_PP[i];
        double expect = correct_return[i];
        double error = actual - expect;
        if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
            printf("* PASSED");
        } else {
            all_passed = 0;
            printf("X FAILED");
        }
        printf(" -> test_vh_calculation at %i",i);
        printf(" -> d1 = %f",d1);
        printf(" -> d2 = %f",d2);
        printf(" -> %f returned",actual);
        printf(" -> %f correct\n",expect);
      }
    } else {
      all_passed = 0;
      printf("X FAILED");
      printf(" -> test_vh_calculation() -- size");
      printf(" -> d1 = %f",d1);
      printf(" -> d2 = %f",d2);
      printf(" -> %f returned",vh_size);
      printf(" -> %f correct\n",correct_size);
    }

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
