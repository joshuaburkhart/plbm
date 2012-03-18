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

double funct(double *d1_d2) {

    double d1,d2;
    d1 = fabs(*(d1_d2));
    d2 = fabs(*(d1_d2+1));

    /*Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);--------------------------------------*/

    //double *A = array_pow(d1,tau1,p,p);
    double Aout[p*p];
    array_pow(Aout,d1,tau1,p,p);
    //double *B = matrx_mlt(2.00,initVh,p,p);
    double Bout[p*p];
    matrx_mlt(Bout,2.00,initVh,p,p);
    //double *C = array_pow(d1,B,p,p);
    double Cout[p*p];
    array_pow(Cout,d1,Bout,p,p);
    //double *D = matrx_sub(1.00,Cout,p,p);
    double Dout[p*p];
    matrx_sub(Dout,1.00,Cout,p,p);
    //double *E = array_mlt(Aout,p,p,Dout);
    double Eout[p*p];
    array_mlt(Eout,Aout,p,p,Dout);
    //double *Vh = array_rdv(Eout,p,p,1.00 - d1*d1);
    double Vhout[p*p];
    array_rdv(Vhout,Eout,p,p,1.00 - d1*d1);

    /*Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);-------------------------------------*/

    //A = array_pow(d2,tau2,q,q);
    double Aout2[q*q];
    array_pow(Aout2,d2,tau2,q,q);
    //B = matrx_mlt(2.00,initVp,q,q);
    double Bout2[q*q];
    matrx_mlt(Bout2,2.00,initVp,q,q);
    //C = array_pow(d2,B,q,q);
    double Cout2[q*q];
    array_pow(Cout2,d2,Bout2,q,q);
    //D = matrx_sub(1.00,Cout2,q,q);
    double Dout2[q*q];
    matrx_sub(Dout2,1.00,Cout2,q,q);
    //E = array_mlt(Aout2,q,q,Dout2);
    double Eout2[q*q];
    array_mlt(Eout2,Aout2,q,q,Dout2);
    //double *Vp = array_rdv(Eout2,q,q,1.00 - d2*d2);
    double Vpout[q*q];
    array_rdv(Vpout,Eout2,q,q,1.00 - d2*d2);

    /*Vh=Vh./det(Vh)^(1/p);-----------------------------------*/

    double dA = matrx_det(Vhout,p);
    double dB = pow(dA,(1.00/((double) p)));
    array_rdv(Vhout,Vhout,p,p,dB);

    /*Vp=Vp./det(Vp)^(1/q);-----------------------------------*/

    dA = matrx_det(Vpout,q);
    dB = pow(dA,(1.00/((double) q)));
    array_rdv(Vpout,Vpout,q,q,dB);

    /*V=kron(Vp,Vh);-----------------------------------*/

    //double *V = kron(Vpout,q,q,Vhout,p,p);
    double Vout[p*p*q*q];
    kron(Vout,Vpout,q,q,Vhout,p,p);

    /*invV=V\eye(n);-----------------------------------*/

    double *invV = matrx_inv(Vout,n);

    /*U=ones(length(X),1);-----------------------------------*/

    double *U = ones(n,1.00);

    /*b=(U'*invV*U)\(U'*invV*X);-----------------------------------*/

    double *A = tran(U,n,1.00);
    double *B = matrx_mlt2(A,1.00,n,invV,n,n);
    double *C = matrx_mlt2(B,1.00,n,X,n,1.00);
    double *D = tran(U,n,1.00);
    double *E = matrx_mlt2(D,1.00,n,invV,n,n);
    double *F = matrx_mlt2(E,1.00,n,U,n,1.00);
    double b = *(C) / *(F); /*should be a 1 x 1 matrix*/

    /*H=X-b;-----------------------------------*/

    double *H = matrx_sub3(X,n,1.00,b);

    /*MSE=(H'*invV*H)/(n-1);-----------------------------------*/

    A = tran(H,n,1.00);
    B = matrx_mlt2(A,1.00,n,invV,n,n);
    C = matrx_mlt2(B,1.00,n,H,n,1.00);
    double MSE = *(C)/(((double) n) - 1.00);

    free(invV);
    //free(Vp);
    //free(V);
    //free(Vh);
    return MSE;
}

void nelmin ( double fn ( double x[] ), int n, double start[], double xmin[],double *ynewlo, double reqmin, double step[], int konvge, int kcount,int *icount, int *numres, int *ifault ) {

    double ccoeff = 0.5;
    double del;
    double dn;
    double dnn;
    double ecoeff = 2.0;
    double eps = 0.001;
    int i;
    int ihi;
    int ilo;
    int j;
    int jcount;
    int l;
    int nn;
    double *p;
    double *p2star;
    double *pbar;
    double *pstar;
    double rcoeff = 1.0;
    double rq;
    double x;
    double *y;
    double y2star;
    double ylo;
    double ystar;
    double z;
    /*TODO: remove input validation*/
    /*  Check the input parameters. */
    if ( reqmin <= 0.0 )
    {
        *ifault = 1;
        return;
    }

    if ( n < 1 )
    {
        *ifault = 1;
        return;
    }

    if ( konvge < 1 )
    {
        *ifault = 1;
        return;
    }
    p = (double *) malloc(n * (n+1) * sizeof(double));
    pstar = (double *) malloc(n * sizeof(double));
    p2star = (double *) malloc(n * sizeof(double));
    pbar = (double *) malloc(n * sizeof(double));
    y = (double *) malloc((n+1) * sizeof(double));
    *icount = 0;
    *numres = 0;
    jcount = konvge;
    dn = ( double ) ( n );
    nn = n + 1;
    dnn = ( double ) ( nn );
    del = 1.0;
    rq = reqmin * dn;
    /*  Initial or restarted loop.*/
    for ( ; ; )
    {
        for ( i = 0; i < n; i++ )
        {
            p[i+n*n] = start[i];
        }
        y[n] = fn ( start );
        *icount = *icount + 1;

        for ( j = 0; j < n; j++ )
        {
            x = start[j];
            start[j] = start[j] + step[j] * del;
            for ( i = 0; i < n; i++ )
            {
                p[i+j*n] = start[i];
            }
            y[j] = fn ( start );
            *icount = *icount + 1;
            start[j] = x;
        }
        /*
        //  The simplex construction is complete.
        //
        //  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
        //  the vertex of the simplex to be replaced.
        */
        ylo = y[0];
        ilo = 0;
        for ( i = 1; i < nn; i++ )
        {
            if ( y[i] < ylo )
            {
                ylo = y[i];
                ilo = i;
            }
        }
        /*
        //  Inner loop.
        */
        for ( ; ; )
        {
            if ( kcount <= *icount )
            {
                break;
            }
            *ynewlo = y[0];
            ihi = 0;

            for ( i = 1; i < nn; i++ )
            {
                if ( *ynewlo < y[i] )
                {
                    *ynewlo = y[i];
                    ihi = i;
                }
            }
            /*
            //  Calculate PBAR, the centroid of the simplex vertices
            //  excepting the vertex with Y value YNEWLO.
            */
            for ( i = 0; i < n; i++ )
            {
                z = 0.0;
                for ( j = 0; j < nn; j++ )
                {
                    z = z + p[i+j*n];
                }
                z = z - p[i+ihi*n];
                pbar[i] = z / dn;
            }
            /*
            //  Reflection through the centroid.
            */
            for ( i = 0; i < n; i++ )
            {
                pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
            }
            ystar = fn ( pstar );
            *icount = *icount + 1;
            /*
            //  Successful reflection, so extension.
            */
            if ( ystar < ylo )
            {
                for ( i = 0; i < n; i++ )
                {
                    p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
                }
                y2star = fn ( p2star );
                *icount = *icount + 1;
                /*
                //  Check extension.
                */
                if ( ystar < y2star )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = pstar[i];
                    }
                    y[ihi] = ystar;
                }
                /*
                //  Retain extension or contraction.
                */
                else
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = p2star[i];
                    }
                    y[ihi] = y2star;
                }
            }
            /*
            //  No extension.
            */
            else
            {
                l = 0;
                for ( i = 0; i < nn; i++ )
                {
                    if ( ystar < y[i] )
                    {
                        l = l + 1;
                    }
                }

                if ( 1 < l )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = pstar[i];
                    }
                    y[ihi] = ystar;
                }
                /*
                //  Contraction on the Y(IHI) side of the centroid.
                */
                else if ( l == 0 )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
                    }
                    y2star = fn ( p2star );
                    *icount = *icount + 1;
                    /*
                    //  Contract the whole simplex.
                    */
                    if ( y[ihi] < y2star )
                    {
                        for ( j = 0; j < nn; j++ )
                        {
                            for ( i = 0; i < n; i++ )
                            {
                                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                                xmin[i] = p[i+j*n];
                            }
                            y[j] = fn ( xmin );
                            *icount = *icount + 1;
                        }
                        ylo = y[0];
                        ilo = 0;

                        for ( i = 1; i < nn; i++ )
                        {
                            if ( y[i] < ylo )
                            {
                                ylo = y[i];
                                ilo = i;
                            }
                        }
                        continue;
                    }
                    /*
                    //  Retain contraction.
                    */
                    else
                    {
                        for ( i = 0; i < n; i++ )
                        {
                            p[i+ihi*n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                }
                /*
                //  Contraction on the reflection side of the centroid.
                */
                else if ( l == 1 )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
                    }
                    y2star = fn ( p2star );
                    *icount = *icount + 1;
                    /*
                    //  Retain reflection?
                    */
                    if ( y2star <= ystar )
                    {
                        for ( i = 0; i < n; i++ )
                        {
                            p[i+ihi*n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                    else
                    {
                        for ( i = 0; i < n; i++ )
                        {
                            p[i+ihi*n] = pstar[i];
                        }
                        y[ihi] = ystar;
                    }
                }
            }
            /*
            //  Check if YLO improved.
            */
            if ( y[ihi] < ylo )
            {
                ylo = y[ihi];
                ilo = ihi;
            }
            jcount = jcount - 1;

            if ( 0 < jcount )
            {
                continue;
            }
            /*
            //  Check to see if minimum reached.
            */
            if ( *icount <= kcount )
            {
                jcount = konvge;

                z = 0.0;
                for ( i = 0; i < nn; i++ )
                {
                    z = z + y[i];
                }
                x = z / dnn;

                z = 0.0;
                for ( i = 0; i < nn; i++ )
                {
                    z = z + pow ( y[i] - x, 2 );
                }

                if ( z <= rq )
                {
                    break;
                }
            }
        }
        /*
        //  Factorial tests to check that YNEWLO is a local minimum.
        */
        for ( i = 0; i < n; i++ )
        {
            xmin[i] = p[i+ilo*n];
        }
        *ynewlo = y[ilo];

        if ( kcount < *icount )
        {
            *ifault = 2;
            break;
        }

        *ifault = 0;

        for ( i = 0; i < n; i++ )
        {
            del = step[i] * eps;
            xmin[i] = xmin[i] + del;
            z = fn ( xmin );
            *icount = *icount + 1;
            if ( z < *ynewlo )
            {
                *ifault = 2;
                break;
            }
            xmin[i] = xmin[i] - del - del;
            z = fn ( xmin );
            *icount = *icount + 1;
            if ( z < *ynewlo )
            {
                *ifault = 2;
                break;
            }
            xmin[i] = xmin[i] + del;
        }

        if ( *ifault == 0 )
        {
            break;
        }
        /*
        //  Restart the procedure.
        */
        for ( i = 0; i < n; i++ )
        {
            start[i] = xmin[i];
        }
        del = eps;
        *numres = *numres + 1;
    }
    free(p);
    free(pstar);
    free(p2star);
    free(pbar);
    free(y);

    return;
}
