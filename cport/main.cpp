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

int main(int argc,char *argv[]) {
    if(argc!=3){
	printf("No arguments supplied\n");
	printf("Usage: a.out .5 .5\n");
	return 1;
    }
    double *d1_d2;
    d1_d2 = (double *) malloc(2 * sizeof(double));
    *(d1_d2) = atof(argv[1]);
    *(d1_d2+1) = atof(argv[2]);

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
