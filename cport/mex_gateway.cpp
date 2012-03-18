#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./lib/lapack.h"
#include "./lib/arraylib.h"
#include "./lib/mex.h"
#include "./lib/minimize.h"

double *initVh;
double *initVp;
double *X;
double *tau1;
double *tau2;
int n;
int p;
int q;

void output(double *matrix,int m,int n);
//void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    /*expect phrs-> 0=initVh, 1=initVp, 2=n, 3=p, 4=q, 5=X, 6=tau1, 7=tau2, 8=d1_d2*/

    /*TODO: n, p, and q may not be necessary as they might already be included in dims of other vars..?*/
    /* http://cnx.org/content/m12348/latest/ */
    /*TODO: n, p, and q are ints -> changing all the doubles may save memory*/

    /* /////////////////// */
    /* convert matlab vals */
    /* /////////////////// */

    /*TODO: these pointers may have to be transposed before asymmetric arrays are properly consumed*/

    n=*mxGetPr(prhs[2]); /* should be 1 x 1 matrix */
    /*also try n=(int) mxGetScalar(prhs[2]); */
    p=*mxGetPr(prhs[3]); /* should be 1 x 1 matrix */
    q=*mxGetPr(prhs[4]); /* should be 1 x 1 matrix */
    initVh=tran(mxGetPr(prhs[0]),p,p);
    initVp=tran(mxGetPr(prhs[1]),q,q);
    X=tran(mxGetPr(prhs[5]),1,n);
    tau1=tran(mxGetPr(prhs[6]),p,p);
    tau2=tran(mxGetPr(prhs[7]),q,q);

    double *d1_d2=mxGetPr(prhs[8]);

    /* //////////////////// */
    /* call nelmin on funct */
    /* //////////////////// */

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

    /* /////////////////////////////// */
    /* return results in matlab format */
    /* /////////////////////////////// */

    /*expects colLen, rowLen, data_type*/
    plhs[0]=mxCreateDoubleMatrix(1,2,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
    double *est = mxGetPr(plhs[0]);
    double *MSE = mxGetPr(plhs[1]);

    est[0]= *(XMIN);
    est[1]= *(XMIN+1);
    MSE[0]= YNEWLO;

    //TODO: free allocated memory ?
    //free(d1_d2);
    free(initVh);
    free(initVp);
    free(X);
    free(tau1);
    free(tau2);
    return;
}

void output(double *matrix,int m,int n) {

    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            mexPrintf("%f ",*(matrix+(i*n+j)));
        }
        mexPrintf("\n");
    }
}
