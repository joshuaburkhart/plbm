#include <stdio.h>
#include <stdlib.h>
#include "jlapack.h"

extern double a_initVh[];
extern double a_initVp[];
extern double a_X[];
extern double a_tau1[];
extern double a_tau2[];
extern double n;
extern double p;
extern double q;

int main(int argc,char *argv[]){
  double d1;
  double d2;
  int i;
  if(argc!=3){
    printf("argument mismatch: you must supply d1 and d2 as doubles...aborting\n");
    return -1;
  }
  if(sizeof(double)==sizeof(argv[1])){
    d1 = atof(argv[1]);
  }else{
    return -2;
  }
  if(sizeof(double)==sizeof(argv[2])){
    d2 = atof(argv[2]);
  }else{
    return -3;
  }
  printf("d1: %f\nd2: %f\n",d1,d2);
  
  output(a_initVh,p,p);
  output(a_initVp,q,q);
  printf("n: %f\n",n);
  printf("p: %f\n",p);
  printf("q: %f\n",q);
  output(a_X,n,1);
  output(a_tau1,p,p);
  output(a_tau2,q,q);

  return 0;
}
