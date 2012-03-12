#include <stdio.h>
#include <stdlib.h>
#include "jlapack.h"

extern double *initVh,initVp,X,tau1,tau2;
extern double n,p,q;

int main(int argc,char *argv[]){
  double d1;
  double d2;
  int i;
  if(argc!=3){
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
  return 0;
}
