#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "arraylib.h"
//#include "mex.h"

double initVh[]= {1.5488,1.2907,0,0.5163,0.7744,0.5163,1.0326,0.2581,0.2581,0.2581,1.0326,0.2581,1.2907,1.5488,0,0.5163,0.7744,0.5163,1.0326,0.2581,0.2581,0.2581,1.0326,0.2581,0,0,1.5488,0,0,0,0,0,0,0,0,0,0.5163,0.5163,0,1.5488,0.5163,0.7744,0.5163,0.2581,0.2581,0.2581,0.5163,0.2581,0.7744,0.7744,0,0.5163,1.5488,0.5163,0.7744,0.2581,0.2581,0.2581,0.7744,0.2581,0.5163,0.5163,0,0.7744,0.5163,1.5488,0.5163,0.2581,0.2581,0.2581,0.5163,0.2581,1.0326,1.0326,0,0.5163,0.7744,0.5163,1.5488,0.2581,0.2581,0.2581,1.2907,0.2581,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,1.5488,0.7744,0.7744,0.2581,0.5163,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,0.7744,1.5488,1.0326,0.2581,0.5163,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,0.7744,1.0326,1.5488,0.2581,0.5163,1.0326,1.0326,0,0.5163,0.7744,0.5163,1.2907,0.2581,0.2581,0.2581,1.5488,0.2581,0.2581,0.2581,0,0.2581,0.2581,0.2581,0.2581,0.5163,0.5163,0.5163,0.2581,1.5488};

double initVp[]= {1.6298,0.6112,0.4074,0.4074,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.6112,1.6298,0.4074,0.4074,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.4074,0.4074,1.6298,0.6112,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.4074,0.4074,0.6112,1.6298,0.2037,0.2037,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,1.6298,0.4074,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,0.4074,1.6298,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.6298,1.2223,1.0186,1.0186,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,1.2223,1.6298,1.0186,1.0186,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,1.0186,1.0186,1.6298,1.2223,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,1.0186,1.0186,1.2223,1.6298,0.8149,0.8149,0.8149,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.8149,0.8149,0.8149,0.8149,1.6298,1.2223,1.0186,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.8149,0.8149,0.8149,0.8149,1.2223,1.6298,1.0186,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.8149,0.8149,0.8149,0.8149,1.0186,1.0186,1.6298,0.6112,0.6112,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,1.6298,0.8149,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.8149,1.6298,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.2037,0.2037,0.4074,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.6298,1.4260,1.2223,1.0186,1.0186,1.0186,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.4260,1.6298,1.2223,1.0186,1.0186,1.0186,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.2223,1.2223,1.6298,1.0186,1.0186,1.0186,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.0186,1.0186,1.0186,1.6298,1.2223,1.2223,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.0186,1.0186,1.0186,1.2223,1.6298,1.2223,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,1.0186,1.0186,1.0186,1.2223,1.2223,1.6298,0.8149,0.8149,0.8149,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.8149,0.8149,0.8149,0.8149,0.8149,0.8149,1.6298,1.2223,1.0186,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.8149,0.8149,0.8149,0.8149,0.8149,0.8149,1.2223,1.6298,1.0186,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.8149,0.8149,0.8149,0.8149,0.8149,0.8149,1.0186,1.0186,1.6298,0.2037,0.2037,0.6112,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,1.6298,0.4074,0.2037,0,0,0,0,0,0,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.2037,0.4074,1.6298,0.2037,0,0,0,0,0,0,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.4074,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.6112,0.2037,0.2037,1.6298};

double X[]= {0,0,0,0,0,0,0,0,0.1885,0.4515,0.2220,0.0662,0,0,0,0,0,0,0.5720,0.5240,0,0,0,0,0,0,0.3816,0.5455,0.6602,0.3370,0,0,0,0,0,0,0.6451,0.5336,0,0,0,0,0,0,0,0,0,0,0,0,0.0658,0.0220,0.0519,0,0.1572,0.2731,0.1463,0.0700,0.1473,0.0445,0.2304,0.1647,0,0,0,0,0,0,0,0,0,0,0.1348,0.1333,0.0917,0.0914,0.1279,0.1228,0.1094,0.0775,0.0830,0.0793,0.2078,0.3226,0.0750,0.0669,0,0.0382,0,0.0383,0.0625,0,0,0,0,0,0,0,0,0.0220,0,0,0,0,0,0,0,0,0,0.0669,0.1249,0.2017,0.1652,0.1224,0,0,0,0,0.0821,0,0.1019,0.1625,0.0911,0.1020,0.0413,0.0956,0,0,0,0,0.0874,0,0.1430,0.1418,0.2610,0.2105,0.1389,0.2819,0.1342,0.3068,0.1372,0.0741,0.1569,0.1360,0,0,0,0.0382,0,0,0,0,0,0,0,0.0402,0,0,0,0,0,0,0,0,0,0.1362,0.1472,0,0.0616,0,0.0882,0.0220,0,0.0951,0,0,0,0,0,0,0.0305,0,0.0825,0.1034,0.1917,0.1462,0,0,0,0,0,0,0.0688,0,0.0510,0.0854,0.0911,0.1852,0.1094,0.1098,0.2535,0.1213,0.2741,0.1146,0.0444,0,0,0,0,0,0,0,0,0.0512,0.1073,0.2325,0,0,0,0.0360,0,0,0,0,0,0,0,0.0548,0.0367,0,0,0.0422,0,0,0,0,0,0.0604,0,0.0827,0.0289,0,0,0.0360,0,0,0,0,0,0,0,0.0548,0.1270,0.1161,0.0856,0.1130,0.2015,0.1348,0.1368,0.2069,0.0830,0.0765,0.2206,0.2203,0,0,0,0,0,0,0,0,0,0,0.0932,0.0922,0.1426,0.1502,0.0349,0.0566,0.0918,0.0383,0,0.1098,0,0.0283,0.0297,0.0400,0.0204,0,0,0,0,0,0.0684,0,0,0,0.0433,0,0.0268,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0551,0,0,0,0,0,0};

double tau1[]= {0,0.5163,3.0977,2.0651,1.5488,2.0651,1.0326,2.5814,2.5814,2.5814,1.0326,2.5814,0.5163,0,3.0977,2.0651,1.5488,2.0651,1.0326,2.5814,2.5814,2.5814,1.0326,2.5814,3.0977,3.0977,0,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,3.0977,2.0651,2.0651,3.0977,0,2.0651,1.5488,2.0651,2.5814,2.5814,2.5814,2.0651,2.5814,1.5488,1.5488,3.0977,2.0651,0,2.0651,1.5488,2.5814,2.5814,2.5814,1.5488,2.5814,2.0651,2.0651,3.0977,1.5488,2.0651,0,2.0651,2.5814,2.5814,2.5814,2.0651,2.5814,1.0326,1.0326,3.0977,2.0651,1.5488,2.0651,0,2.5814,2.5814,2.5814,0.5163,2.5814,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,0,1.5488,1.5488,2.5814,2.0651,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,1.5488,0,1.0326,2.5814,2.0651,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,1.5488,1.0326,0,2.5814,2.0651,1.0326,1.0326,3.0977,2.0651,1.5488,2.0651,0.5163,2.5814,2.5814,2.5814,0,2.5814,2.5814,2.5814,3.0977,2.5814,2.5814,2.5814,2.5814,2.0651,2.0651,2.0651,2.5814,0};

double tau2[]= {0,2.0372,2.4447,2.4447,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.0372,0,2.4447,2.4447,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,0,2.0372,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.0372,0,2.8521,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,0,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,2.4447,0,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,0,0.8149,1.2223,1.2223,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,0.8149,0,1.2223,1.2223,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.2223,1.2223,0,0.8149,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.2223,1.2223,0.8149,0,1.6298,1.6298,1.6298,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.6298,1.6298,1.6298,1.6298,0,0.8149,1.2223,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.6298,1.6298,1.6298,1.6298,0.8149,0,1.2223,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,1.6298,1.6298,1.6298,1.6298,1.2223,1.2223,0,2.0372,2.0372,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,0,1.6298,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,1.6298,0,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.8521,2.8521,2.4447,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,0,0.4074,0.8149,1.2223,1.2223,1.2223,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,0.4074,0,0.8149,1.2223,1.2223,1.2223,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,0.8149,0.8149,0,1.2223,1.2223,1.2223,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.2223,1.2223,1.2223,0,0.8149,0.8149,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.2223,1.2223,1.2223,0.8149,0,0.8149,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.2223,1.2223,1.2223,0.8149,0.8149,0,1.6298,1.6298,1.6298,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.6298,1.6298,1.6298,1.6298,1.6298,1.6298,0,0.8149,1.2223,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.6298,1.6298,1.6298,1.6298,1.6298,1.6298,0.8149,0,1.2223,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,1.6298,1.6298,1.6298,1.6298,1.6298,1.6298,1.2223,1.2223,0,2.8521,2.8521,2.0372,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,0,2.4447,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.8521,2.4447,0,2.8521,3.2595,3.2595,3.2595,3.2595,3.2595,3.2595,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.4447,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.0372,2.8521,2.8521,0};
double n=324; //this is p * q... ?
double p=12;
double q=27;

double funct(double *d1_d2);
void nelmin ( double fn ( double x[] ), int n, double start[], double xmin[],double *ynewlo, double reqmin, double step[], int konvge, int kcount,int *icount, int *numres, int *ifault );

/*
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){

    //should grab all the globs from matlab, set vals, call nelmin on funct, and return results in matlab format

    double STEP[2];
    if(*(d1_d2)==0) {
        STEP[0]=0.00025;
    } else {
        STEP[0]=0.95 * *(d1_d2);
    }
    if(*(d1_d2+1)==0) {
        STEP[1]=0.00025;
    } else {
        STEP[1]=0.95 * *(d1_d2+1);
    }
    double XMIN[2]; //coordinates of minimum value
    double YNEWLO; //minimum value
    double REQMIN = 0.0001; //termination variance limit
    int KONVGE = 10; //frequency of convergence tests
    int KCOUNT = 10000; //max number of iterations
    int ICOUNT; //number of evaluations
    int NUMRES; //number of restarts
    int IFAULT; //error indicator
    nelmin(funct,2,d1_d2,XMIN,&YNEWLO,REQMIN,STEP,KONVGE,KCOUNT,&ICOUNT,&NUMRES,&IFAULT);

    printf("minimization coordinates: %f %f\n",*(XMIN),*(XMIN+1));
    printf("minimum value: %f\n",YNEWLO);

    free(d1_d2);
  return;
}
*/

int main(int argc,char *argv[]) {
    double *d1_d2;
    int i;
    d1_d2 = (double *) malloc(2 * sizeof(double));
    *(d1_d2) = atof(argv[1]);
    *(d1_d2+1) = atof(argv[2]);

    double STEP[2];
    if(*(d1_d2)==0) {
        STEP[0]=0.00025;
    } else {
        STEP[0]=0.95 * *(d1_d2);
    }
    if(*(d1_d2+1)==0) {
        STEP[1]=0.00025;
    } else {
        STEP[1]=0.95 * *(d1_d2+1);
    }
    double XMIN[2]; //coordinates of minimum value
    double YNEWLO; //minimum value
    double REQMIN = 0.0001; //termination variance limit
    int KONVGE = 10; //frequency of convergence tests
    int KCOUNT = 10000; //max number of iterations
    int ICOUNT; //number of evaluations
    int NUMRES; //number of restarts
    int IFAULT; //error indicator
    nelmin(funct,2,d1_d2,XMIN,&YNEWLO,REQMIN,STEP,KONVGE,KCOUNT,&ICOUNT,&NUMRES,&IFAULT);

    printf("minimization coordinates: %f %f\n",*(XMIN),*(XMIN+1));
    printf("minimum value: %f\n",YNEWLO);

    free(d1_d2);
    return 0;
}

double funct(double *d1_d2) {

    double d1,d2;
    d1 = fabs(*(d1_d2));
    d2 = fabs(*(d1_d2+1));

    //Vh=(d1.^tau1).*(1-d1.^(2*initVh))./(1-d1^2);--------------------------------------Vh

    double *A = array_pow(d1,tau1,p,p);
    double *B = matrx_mlt(2,initVh,p,p);
    double *C = array_pow(d1,B,p,p);
    double *D = matrx_sub(1,C,p,p);
    double *E = array_mlt(A,p,p,D);
    double d1sq = d1*d1;
    double omd1 = 1 - d1sq;
    double *Vh = array_rdv(E,p,p,omd1);

    //Vp=(d2.^tau2).*(1-d2.^(2*initVp))./(1-d2^2);-------------------------------------Vp

    A = array_pow(d2,tau2,q,q);
    B = matrx_mlt(2,initVp,q,q);
    C = array_pow(d2,B,q,q);
    D = matrx_sub(1,C,q,q);
    E = array_mlt(A,q,q,D);
    double d2sq = d2*d2;
    double omd2 = 1 - d2sq;
    double *Vp = array_rdv(E,q,q,omd2);

    //Vh=Vh./det(Vh)^(1/p);-----------------------------------Vh

    double dtm = matrx_det(Vh,p);
    double dtm21op = pow(dtm,(1/p));
    Vh = array_rdv(Vh,p,p,dtm21op);

    //Vp=Vp./det(Vp)^(1/q);-----------------------------------Vp

    dtm = matrx_det(Vp,q);
    double dtm21oq = pow(dtm,(1/q));
    Vp = array_rdv(Vp,q,q,dtm21oq);

    //V=kron(Vp,Vh);-----------------------------------V

    double *V = kron(Vp,q,q,Vh,p,p);

    //invV=V\eye(n);-----------------------------------invV

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
        printf("d1: %f\n",d1);
        printf("d2: %f\n",d2);
        info=0;
    }
    dgetri_(&N,A,&lda,ipiv,work,&lwork,&info);
    if(info!=0) {
        printf("dgetri returns info code %i\n",info);
        printf("d1: %f\n",d1);
        printf("d2: %f\n",d2);
    }

    invV = tran(A,n,n); //column major -> row major

    //U=ones(length(X),1);-----------------------------------U

    double *U = ones(n,1);

    //b=(U'*invV*U)\(U'*invV*X);-----------------------------------b

    A = tran(U,n,1);
    B = matrx_mlt2(A,1,n,invV,n,n);
    C = matrx_mlt2(B,1,n,U,n,1);
    D = matrx_mlt2(B,1,n,X,n,1);

    double c = *(C); //should be a 1 x 1 matrix
    double d = *(D); //should be a 1 x 1 matrix
    double b = d/c;

    //H=X-b;-----------------------------------H

    double *H = matrx_sub3(X,n,1,b);

    //MSE=(H'*invV*H)/(n-1);-----------------------------------MSE

    A = tran(H,n,1);
    B = matrx_mlt2(A,1,n,invV,n,n);
    C = matrx_mlt2(B,1,n,H,n,1);
    c = *(C); //should be a 1 x 1 matrix
    double MSE = c/(n -1);

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
    //
    //  Check the input parameters.
    //
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
    //
    //  Initial or restarted loop.
    //
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
        //
        //  The simplex construction is complete.
        //
        //  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
        //  the vertex of the simplex to be replaced.
        //
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
        //
        //  Inner loop.
        //
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
            //
            //  Calculate PBAR, the centroid of the simplex vertices
            //  excepting the vertex with Y value YNEWLO.
            //
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
            //
            //  Reflection through the centroid.
            //
            for ( i = 0; i < n; i++ )
            {
                pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
            }
            ystar = fn ( pstar );
            *icount = *icount + 1;
            //
            //  Successful reflection, so extension.
            //
            if ( ystar < ylo )
            {
                for ( i = 0; i < n; i++ )
                {
                    p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
                }
                y2star = fn ( p2star );
                *icount = *icount + 1;
                //
                //  Check extension.
                //
                if ( ystar < y2star )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = pstar[i];
                    }
                    y[ihi] = ystar;
                }
                //
                //  Retain extension or contraction.
                //
                else
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p[i+ihi*n] = p2star[i];
                    }
                    y[ihi] = y2star;
                }
            }
            //
            //  No extension.
            //
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
                //
                //  Contraction on the Y(IHI) side of the centroid.
                //
                else if ( l == 0 )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
                    }
                    y2star = fn ( p2star );
                    *icount = *icount + 1;
                    //
                    //  Contract the whole simplex.
                    //
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
                    //
                    //  Retain contraction.
                    //
                    else
                    {
                        for ( i = 0; i < n; i++ )
                        {
                            p[i+ihi*n] = p2star[i];
                        }
                        y[ihi] = y2star;
                    }
                }
                //
                //  Contraction on the reflection side of the centroid.
                //
                else if ( l == 1 )
                {
                    for ( i = 0; i < n; i++ )
                    {
                        p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
                    }
                    y2star = fn ( p2star );
                    *icount = *icount + 1;
                    //
                    //  Retain reflection?
                    //
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
            //
            //  Check if YLO improved.
            //
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
            //
            //  Check to see if minimum reached.
            //
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
        //
        //  Factorial tests to check that YNEWLO is a local minimum.
        //
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
        //
        //  Restart the procedure.
        //
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