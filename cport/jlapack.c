#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "jlapack.h"
//#include <clapack.h>

double a_initVh[]= {1.5488,1.2907,0,0.5163,0.7744,0.5163,1.0326,0.2581,0.2581,0.2581,1.0326,0.2581,1.2907,1.5488,0,0.5163,0.7744,0.5163,1.0326,0.2581,0.2581,0.2581,1.0326,0.2581,0,0,1.5488,0,0,0,0,0,0,0,0,0,0.5163,0.5163,0,1.5488,0.5163,0.7744,0.5163,0.2581,0.2581,0.2581,0.5163,0.2581,0.7744,0.7744,0,0.5163,1.5488,0.5163,0.7744,0.2581,0.2581,0.2581,0.7744,0.2581,0.5163,0.5163,0,0.7744,0.5163,1.5488,0.5163,0.2581,0.2581,0.2581,0.5163,0.2581,1.0326,1.0326,0,0.5163,0.7744,0.5163,1.5488,0.2581,0.2581,0.2581,1.2907,0.2581,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,1.5488,0.7744,0.7744,0.2581,0.5163,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,0.7744,1.5488,1.0326,0.2581,0.5163,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,0.7744,1.0326,1.5488,0.2581,0.5163,1.0326,1.0326,0,0.5163,0.7744,0.5163,1.2907,0.2581,0.2581,0.2581,1.5488,0.2581,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,0.5163,0.5163,0.5163,0.2581,1.5488};

double a_initVp[]= {1.6298,0.6112,0.4074,0.4074,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.6112,1.6298,0.4074,0.4074,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.4074,0.4074,1.6298,0.6112,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.4074,0.4074,0.6112,1.6298,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,1.6298,0.4074,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,0.4074,1.6298,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.6298,1.2223,1.0186,1.0186,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,1.2223,1.6298,1.0186,1.0186,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,1.0186,1.0186,1.6298,1.2223,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,1.0186,1.0186,1.2223,1.6298,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.8149,0.8149,0.8149,0.8149,1.6298,1.2223,1.0186,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.8149,0.8149,0.8149,0.8149,1.2223,1.6298,1.0186,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.8149,0.8149,0.8149,0.8149,1.0186,1.0186,1.6298,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,1.6298,0.8149,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.8149,1.6298,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.6298,1.4260,1.2223,1.0186,1.0186,1.0186,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.4260,1.6298,1.2223,1.0186,1.0186,1.0186,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.2223,1.2223,1.6298,1.0186,1.0186,1.0186,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.0186,1.0186,1.0186,1.6298,1.2223,1.2223,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.0186,1.0186,1.0186,1.2223,1.6298,1.2223,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.0186,1.0186,1.0186,1.2223,1.2223,1.6298,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.8149,0.8149,0.8149,0.8149,0.8149,0.8149,1.6298,1.2223,1.0186,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.8149,0.8149,0.8149,0.8149,0.8149,0.8149,1.2223,1.6298,1.0186,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.8149,0.8149,0.8149,0.8149,0.8149,0.8149,1.0186,1.0186,1.6298,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,1.6298,0.4074,0.2037,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.4074,1.6298,0.2037,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.2037,0.2037,1.6298};

double n=324; //this is p * q... ?
double p=12;
double q=27;

double a_X[]= {0,0,0,0,0,0,0,0,0.1885,0.4515,0.2220,0.0662,0,0,0,0,0,0,0.5720,0.5240,0,0,0,0,0,0,0.3816,0.5455,0.6602,0.3370,0,0,0,0,0,0,0.6451,0.5336,0,0,0,0,0,0,0,0,0,0,0,0,0.0658,0.0220,0.0519,0,0.1572,0.2731,0.1463,0.0700,0.1473,0.0445,0.2304,0.1647,0,0,0,0,0,0,0,0,0,0,0.1348,0.1333,0.0917,0.0914,0.1279,0.1228,0.1094,0.0775,0.0830,0.0793,0.2078,0.3226,0.0750,0.0669,0,0.0382,0,0.0383,0.0625,0,0,0,0,0,0,0,0,0.0220,0,0,0,0,0,0,0,0,0,0.0669,0.1249,0.2017,0.1652,0.1224,0,0,0,0,0.0821,0,0.1019,0.1625,0.0911,0.1020,0.0413,0.0956,0,0,0,0,0.0874,0,0.1430,0.1418,0.2610,0.2105,0.1389,0.2819,0.1342,0.3068,0.1372,0.0741,0.1569,0.1360,0,0,0,0.0382,0,0,0,0,0,0,0,0.0402,0,0,0,0,0,0,0,0,0,0.1362,0.1472,0,0.0616,0,0.0882,0.0220,0,0.0951,0,0,0,0,0,0,0.0305,0,0.0825,0.1034,0.1917,0.1462,0,0,0,0,0,0,0.0688,0,0.0510,0.0854,0.0911,0.1852,0.1094,0.1098,0.2535,0.1213,0.2741,0.1146,0.0444,0,0,0,0,0,0,0,0,0.0512,0.1073,0.2325,0,0,0,0.0360,0,0,0,0,0,0,0,0.0548,0.0367,0,0,0.0422,0,0,0,0,0,0.0604,0,0.0827,0.0289,0,0,0.0360,0,0,0,0,0,0,0,0.0548,0.1270,0.1161,0.0856,0.1130,0.2015,0.1348,0.1368,0.2069,0.0830,0.0765,0.2206,0.2203,0,0,0,0,0,0,0,0,0,0,0.0932,0.0922,0.1426,0.1502,0.0349,0.0566,0.0918,0.0383,0,0.1098,0,0.0283,0.0297,0.0400,0.0204,0,0,0,0,0,0.0684,0,0,0,0.0433,0,0.0268,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0551,0,0,0,0,0,0};

double a_tau1[]= {0,0.5163,3.0977,2.0651,1.5488,2.0651,1.0326,2.5814,2.5814,2.5814,1.0326,2.5814,0.5163,0,3.0977,2.0651,1.5488,2.0651,1.0326,2.5814,2.5814,2.5814,1.0326,2.5814,3.0977,3.0977,0,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,2.0651,2.0651,3.0977,0,2.0651,1.5488,2.0651,2.5814,2.5814,2.5814,2.0651,2.5814,1.5488,1.5488,3.0977,2.0651,0,2.0651,1.5488,2.5814,2.5814,2.5814,1.5488,2.5814,2.0651,2.0651,3.0977,1.5488,2.0651,0,2.0651,2.5814,2.5814,2.5814,2.0651,2.5814,1.0326,1.0326,3.0977,2.0651,1.5488,2.0651,0,2.5814,2.5814,2.5814,0.5163,2.5814,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,0,1.5488,1.5488,2.5814,2.0651,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,1.5488,0,1.0326,2.5814,2.0651,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,1.5488,1.0326,0,2.5814,2.0651,1.0326,1.0326,3.0977,2.0651,1.5488,2.0651,0.5163,2.5814,2.5814,2.5814,0,2.5814,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,2.0651,2.0651,2.0651,2.5814,0};

double a_tau2[]= {0,2.0372,2.4447,2.4447,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.0372,0,2.4447,2.4447,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,0,2.0372,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.0372,0,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,0,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,2.4447,0,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,0,0.8149,1.2223,1.2223,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,0.8149,0,1.2223,1.2223,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.2223,1.2223,0,0.8149,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.2223,1.2223,0.8149,0,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.6298,1.6298,1.6298,1.6298,0,0.8149,1.2223,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.6298,1.6298,1.6298,1.6298,0.8149,0,1.2223,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.6298,1.6298,1.6298,1.6298,1.2223,1.2223,0,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,0,1.6298,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,1.6298,0,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,0,0.4074,0.8149,1.2223,1.2223,1.2223,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,0.4074,0,0.8149,1.2223,1.2223,1.2223,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,0.8149,0.8149,0,1.2223,1.2223,1.2223,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.2223,1.2223,1.2223,0,0.8149,0.8149,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.2223,1.2223,1.2223,0.8149,0,0.8149,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.2223,1.2223,1.2223,0.8149,0.8149,0,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.6298,1.6298,1.6298,1.6298,1.6298,1.6298,0,0.8149,1.2223,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.6298,1.6298,1.6298,1.6298,1.6298,1.6298,0.8149,0,1.2223,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.6298,1.6298,1.6298,1.6298,1.6298,1.6298,1.2223,1.2223,0,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,0,2.4447,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.4447,0,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.8521,2.8521,0};

/* main is for testing purposes
int main(void){

  double *result;

  printf("---------------identity\n");

  result = eye(L);
  output(result,L,L);

  printf("---------------kronecker\n");

  double A[4]={1,2,3,4};
  double B[4]={0,5,6,7};
  result=kron(A,2,2,B,2,2);
  output(result,4,4);

  printf("---------------transpose\n");

  double C[6]={1,2,3,4,5,6};
  result=tran(C,2,3);
  output(result,3,2);

  printf("---------------ones\n");

  result=ones(5,7);
  output(result,5,7);

  printf("---------------determinant\n");

  double d=matrx_det(A,2);
  printf("%f\n",d);

  printf("---------------array_pow\n");

  result=array_pow(5.00,A,2,2);
  output(result,2,2);

  printf("---------------array_mlt\n");

  result=array_mlt(A,2,2,B);
  output(result,2,2);

  printf("---------------matrx_mlt\n");

  result=matrx_mlt(5.00,A,2,2);
  output(result,2,2);

  printf("---------------matrx_mlt2\n");

  result=matrx_mlt2(A,2,2,B,2,2);
  output(result,2,2);

  printf("---------------matrx_sub\n");

  result=matrx_sub(5.00,A,2,2);
  output(result,2,2);

  printf("---------------matrx_sub2\n");

  result=matrx_sub2(A,2,2,B);
  output(result,2,2);

  printf("---------------array_rdv\n");

  result=array_rdv(A,2,2,5.00);
  output(result,2,2);

  free(result);
  return 0;
}
*/

double* array_rdv(double *A,int m,int n,double d) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=*(A+(i*n+j)) / d;
        }
    }
    return diff;
}

double* matrx_sub3(double *A,int m,int n,double d) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=*(A+(i*n+j)) - d;
        }
    }
    return diff;
}


double* matrx_sub2(double *A,int m,int n,double *B) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=*(A+(i*n+j)) - *(B+(i*n+j));
        }
    }
    return diff;
}

double* matrx_sub(double d,double *A,int m,int n) {

    double *diff;
    diff=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(diff+(i*n+j))=d - *(A+(i*n+j));
        }
    }
    return diff;
}

double* matrx_mlt2(double *A,int ma,int na,double *B,int mb,int nb) {

    double *pdct;
    pdct=(double *) malloc(ma*nb*sizeof(double));
    int i;
    int j;
    int k;
    for(i=0; i<ma; i++) {
        for(j=0; j<nb; j++) {
            double sum=0;
            for(k=0; k<na; k++) {
                sum+=*(A+(i*na+k)) * *(B+(j+k*nb));
            }
            *(pdct+(i*nb+j))=sum;
        }
    }
    return pdct;
}

double* matrx_mlt(double d,double *A,int m,int n) {

    double *pdct;
    pdct=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(pdct+(i*n+j))=*(A+(i*n+j)) * d;
        }
    }
    return pdct;
}

double* array_mlt(double *A,int m,int n,double *B) {

    double *pdct;
    pdct=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(pdct+(i*n+j))=*(A+(i*n+j)) * *(B+(i*n+j));
        }
    }
    return pdct;
}

double* array_pow(double d,double *A,int m,int n) {

    double *pdct;
    pdct=(double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            if(d>0) {
                *(pdct+(i*n+j))=pow(d,*(A+(i*n+j)));
            }
            else {
                *(pdct+(i*n+j))=-1 * pow(abs(d),*(A+(i*n+j)));
            }
        }
    }
    return pdct;
}

double matrx_det(double *A,int n) {
    A = tran(A,n,n); //row major -> column major
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
    A = tran(A,n,n);
    double diag=1;
    int i;
    for(i = 0; i < n; i++) {
        double multiplier =  *(A + (i * n + i));
        diag *= multiplier;
	printf("diag: %f\n",diag);
    }
    for(i = 0; i < n; i++) {
        *(A + (i *  n + i)) = 1;
        int j;
        for(j = i+1; j < n; j++) {
            *(A + ( i *  n + j)) = 0;
        }
    }
    double dtm=det_l(A,n);
    printf("dtm from internal: %f\n",dtm);
    return(dtm * diag);
}


//reference: http://cboard.cprogramming.com/cplusplus-programming/30001-determinant-calculation.html
static double det_l(double *A,int n) {
    int i, j, k;
    double **m;
    double det = 1;
    m = (double **) malloc(n*sizeof(double *));
    for ( i = 0; i < n; i++ ) {
        m[i] = (double *) malloc(n*sizeof(double));
    }
    for ( i = 0; i < n; i++ ) {
        for ( j = 0; j < n; j++ ) {
            m[i][j] = *(A+(i*n+j));
        }
    }
    for ( k = 0; k < n; k++ ) {
        if ( m[k][k] == 0 ) {
            int ok = 0;
            for ( j = k; j < n; j++ ) {
                if (m[j][k] != 0 ) {
                    ok = 1;
                }
            }
            if (ok==0) {
                return 0;
            }
            for ( i = k; i < n; i++ ) {
                double tmp= m[i][j];
                m[i][j]=m[i][k];
                m[i][k]=tmp;
            }
            det = -det;
        }
        det *= m[k][k];
        if ( k + 1 < n ) {
            for ( i = k + 1; i < n; i++ ) {
                for ( j = k + 1; j < n; j++ ) {
                    m[i][j] = m[i][j] - m[i][k] * m[k][j] / m[k][k];
                }
            }
        }
    }
    for ( i = 0; i < n; i++ ) {
        free(m[i]);
    }
    free(m);
    return det;
}

double* ones(int m,int n) {

    double *o;
    o = (double *) malloc(m*n*sizeof(double));
    int i;
    for(i=0; i<m*n; i++) {
        *(o+i)=1;
    }
    return o;
}

double* eye(int n) {

    double *iden;
    iden = (double *) malloc(n * n * sizeof(double));
    memset(iden,0,n*n*sizeof(double));
    int i;
    for(i=0; i<n; i++) {
        *(iden+(i*n+i))=1;
    }
    return iden;
}

double* kron(double *A,int ma,int na,double *B,int mb,int nb) {

    double* k;
    k = (double *) malloc(ma*mb*na*nb*sizeof(double));
    int i;
    int j;
    int min_m = ma < mb ? ma : mb;
    int min_n = na < nb ? na : nb;
    for(i=0; i<ma*mb; i++) {
        for(j=0; j<na*nb; j++) {
            int a_row = i/min_m;
            int a_col = j/min_n;
            int b_row = i%min_m;
            int b_col = j%min_n;
            //printf("A[%i][%i],B[%i][%i] ",a_row,a_col,b_row,b_col);
            double val=*(A+(a_row*na+a_col)) * *(B+(b_row*nb+b_col));
            *(k+((i*na*nb)+j))=val;
        }
        //printf("\n");
    }
    return k;
}

double* tran(double *A,int m,int n) {

    double* t;
    t = (double *) malloc(m*n*sizeof(double));
    int i;
    int j;
    for(i=0; i<m; i++) {
        for(j=0; j<n; j++) {
            *(t+(i+j*m))=*(A+(i*n+j));
        }
    }
    return t;
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
