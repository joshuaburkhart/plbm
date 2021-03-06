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
void test_vp_calculation(double d1, double d2, double correct_return[], int correct_size, double epsilon);
void test_vh_rdv_calculation(double io_vh[], double correct_return[], int correct_size, double epsilon);
void test_matrx_det(double A[],double expect,int size, double epsilon);
void test_matrx_mlt(double io[], double correct_return[], double d, double m, double n, double epsilon);
void test_array_pow(double io[],double correct_return[], double d,double m,double n, double epsilon);
void test_matrx_sub(double io[],double correct_return[], double d,double m,double n,double epsilon);
void test_array_amlt(double io[],double correct_return[], double b[],double m,double n,double epsilon);
void test_array_rdv(double io[],double correct_return[],double d,double m,double n,double epsilon);
void test_kron(double io[],double correct_return[],double b[],double ma,double na,double mb,double nb,double epsilon);
void test_tran(double io[],double correct_return[],double m,double n,double epsilon);
void test_matrx_inv_(double io[],double correct_return[],int n,double epsilon);
void test_ones(double io[],double correct_return[],int m,int n,double epsilon);
void test_matrx_mlt2(double io[],double io1[],int m,int n,double io2[],int ma,int na,double correct_return[],double epsilon);
void test_matrx_sub3(double io[],double correct_return[], double d,double m,double n,double epsilon);


int main(int argc,char *argv[]) {

	printf("\n\ttesting...\n\n");
	double epsilon=0.00000001;
	double big_epsilon=0.0005;

	//test arraylib functions

	//matrx_mlt

	double mlt_test_input_1[]={1,2,3,4};
	double mlt_test_input_d_1=3;
	double mlt_expect_1[]={3,6,9,12};
	int mlt_length_on_side_1 = sqrt(sizeof(mlt_test_input_1) / sizeof(double));
	test_matrx_mlt(mlt_test_input_1, mlt_expect_1, mlt_test_input_d_1, mlt_length_on_side_1, mlt_length_on_side_1, epsilon);

	double mlt_test_input_2[]={-0.1,2,-0.3,4};
	double mlt_test_input_d_2=-5;
	double mlt_expect_2[]={0.5,-10.0,1.5,-20.0};
	int mlt_length_on_side_2 = sqrt(sizeof(mlt_test_input_2) / sizeof(double));
	test_matrx_mlt(mlt_test_input_2, mlt_expect_2, mlt_test_input_d_2, mlt_length_on_side_2, mlt_length_on_side_2, epsilon);

	double mlt_test_input_3[]={-0.1,2,-0.3,4};
	double mlt_test_input_d_3=0;
	double mlt_expect_3[]={0,0,0,0};
	int mlt_length_on_side_3 = sqrt(sizeof(mlt_test_input_3) / sizeof(double));
	test_matrx_mlt(mlt_test_input_3, mlt_expect_3, mlt_test_input_d_3, mlt_length_on_side_3, mlt_length_on_side_3, epsilon);

	double mlt_test_input_4[]={0,0,0,0};
	double mlt_test_input_d_4=5;
	double mlt_expect_4[]={0,0,0,0};
	int mlt_length_on_side_4 = sqrt(sizeof(mlt_test_input_4) / sizeof(double));
	test_matrx_mlt(mlt_test_input_4, mlt_expect_4, mlt_test_input_d_4, mlt_length_on_side_4, mlt_length_on_side_4, epsilon);

	double mlt_test_input_5[]={0,2,4,6};
	double mlt_test_input_d_5=0.5;
	double mlt_expect_5[]={0,1,2,3};
	int mlt_length_on_side_5 = sqrt(sizeof(mlt_test_input_5) / sizeof(double));
	test_matrx_mlt(mlt_test_input_5, mlt_expect_5, mlt_test_input_d_5, mlt_length_on_side_5, mlt_length_on_side_5, epsilon);

	//array_pow

	double pow_test_input_1[]={-0.1,2,-0.3,4};
	double pow_test_input_d_1=2;
	double pow_expect_1[]={0.9330,4,0.8123,16};
	int pow_length_on_side_1= sqrt(sizeof(pow_test_input_1) / sizeof(double));
	test_array_pow(pow_test_input_1,pow_expect_1,pow_test_input_d_1,pow_length_on_side_1,pow_length_on_side_1,epsilon);

	double pow_test_input_2[]={-0.1,2,-0.3,4};
	double pow_test_input_d_2=3;
	double pow_expect_2[]={0.8960,9,0.7192,81};
	int pow_length_on_side_2= sqrt(sizeof(pow_test_input_2) / sizeof(double));
	test_array_pow(pow_test_input_2,pow_expect_2,pow_test_input_d_2,pow_length_on_side_2,pow_length_on_side_2,epsilon);

	double pow_test_input_3[]={-1,0,1,2};
	double pow_test_input_d_3=2;
	double pow_expect_3[]={0.5,1.0,2.0,4.0};
	int pow_length_on_side_3= sqrt(sizeof(pow_test_input_3) / sizeof(double));
	test_array_pow(pow_test_input_3,pow_expect_3,pow_test_input_d_3,pow_length_on_side_3,pow_length_on_side_3,epsilon);

	double pow_test_input_4[]={-1,0,1,2};
	double pow_test_input_d_4=3;
	double pow_expect_4[]={0.3333,1.0,3.0,9.0};
	int pow_length_on_side_4= sqrt(sizeof(pow_test_input_4) / sizeof(double));
	test_array_pow(pow_test_input_4,pow_expect_4,pow_test_input_d_4,pow_length_on_side_4,pow_length_on_side_4,epsilon);

	double pow_test_input_5[]={-1,0,1,2};
	double pow_test_input_d_5=1;
	double pow_expect_5[]={1,1,1,1};
	int pow_length_on_side_5= sqrt(sizeof(pow_test_input_5) / sizeof(double));
	test_array_pow(pow_test_input_5,pow_expect_5,pow_test_input_d_5,pow_length_on_side_5,pow_length_on_side_5,epsilon);

	//matrx_sub

	double sub_test_input_1[]={1,2,3,4};
	double sub_test_input_d_1=5;
	double sub_expect_1[]={4,3,2,1};
	int sub_length_on_side_1= sqrt(sizeof(sub_test_input_1) / sizeof(double));
	test_matrx_sub(sub_test_input_1,sub_expect_1,sub_test_input_d_1,sub_length_on_side_1,sub_length_on_side_1,epsilon);

	//array_mlt

	double amlt_test_input_1[]={1,2,3,4};
	double amlt_test_input_b_1[]={5,5,5,5};
	double amlt_expect_1[]={5,10,15,20};
	int amlt_length_on_side_1= sqrt(sizeof(amlt_test_input_1) / sizeof(double));
	test_array_amlt(amlt_test_input_1,amlt_expect_1,amlt_test_input_b_1,amlt_length_on_side_1,amlt_length_on_side_1,epsilon);

	//array_rdv

	double ardv_test_input_1[]={1,2,3,4};
	double ardv_test_input_d_1=4;
	double ardv_expect_1[]={0.25,0.5,0.75,1.0};
	int ardv_length_on_side_1= sqrt(sizeof(ardv_test_input_1) / sizeof(double));
	test_array_rdv(ardv_test_input_1,ardv_expect_1,ardv_test_input_d_1,ardv_length_on_side_1,ardv_length_on_side_1,epsilon);

	//kron

	double kron_test_input_1[]={1,2,3,4};
	double kron_test_input_b_1[]={3};
	double kron_expect_1[]={3,6,9,12};
	int kron_length_on_side_a_1= sqrt(sizeof(kron_test_input_1) / sizeof(double));
	int kron_length_on_side_b_1= sqrt(sizeof(kron_test_input_b_1) / sizeof(double));
	test_kron(kron_test_input_1,kron_expect_1,kron_test_input_b_1,kron_length_on_side_a_1,kron_length_on_side_a_1,kron_length_on_side_b_1,kron_length_on_side_b_1,epsilon);

	//tran

	double tran_test_input_1[]={2,3,4,5,6,7,8,9,0};
	double tran_expect_1[]={2,5,8,3,6,9,4,7,0};
	int tran_length_on_side_1= sqrt(sizeof(tran_test_input_1) / sizeof(double));
	test_tran(tran_test_input_1,tran_expect_1,tran_length_on_side_1,tran_length_on_side_1,epsilon);

        double tran_test_input_2[]={1,2,3,4,5,6,7,8};
	double tran_expect_2[]={1,2,3,4,5,6,7,8};
	int tran_num_rows_2=8;
	int tran_num_cols_2=1;
	test_tran(tran_test_input_2,tran_expect_2,tran_num_rows_2,tran_num_cols_2,epsilon);

	//matrx_inv

	double matrx_inv_test_input_1[]={1,2,2,3};
	double matrx_inv_expect_1[]={-3,2,2,-1};
	int matrx_inv_length_on_side_1= sqrt(sizeof(matrx_inv_test_input_1) / sizeof(double));
	test_matrx_inv_(matrx_inv_test_input_1,matrx_inv_expect_1,matrx_inv_length_on_side_1,epsilon);

	double matrx_inv_test_input_2[]={-1.0000,0.5000,0,-0.4000,0.5000,1.0000,1.0000,1.0000,0,1.0000,1.0000,1.0000,-0.4000,1.0000,1.0000,9.0000};
	double matrx_inv_expect_2[]={0.0000,2.0000,-2.0000,-0.0000,2.0000,4.0800,-4.1800,0.1000,-2.0000,-4.1800,5.4050,-0.2250,-0.0000,0.1000,-0.2250,0.1250};
	int matrx_inv_length_on_side_2= sqrt(sizeof(matrx_inv_test_input_2) / sizeof(double));
	test_matrx_inv_(matrx_inv_test_input_2,matrx_inv_expect_2,matrx_inv_length_on_side_2,epsilon);

	//ones

	double ones_test_input_1[]={0,0,0,0,0,0,0,0,0,0};
	double ones_expect_1[]={1,1,1,1,1,1,1,1,1};
	int ones_length_on_side_1= sqrt(sizeof(ones_test_input_1) / sizeof(double));
	test_ones(ones_test_input_1,ones_expect_1,ones_length_on_side_1,ones_length_on_side_1,epsilon);  

	//matrx_mlt2

	double matrx_mlt2_test_input_1[]={2,3,4,5,6,7,8,9,0};
	double matrx_mlt2_test_input_b_1[]={1,2,3};
	double matrx_mlt2_expect_1[]={20,38,26};
	double matrx_mlt2_out_1[]={0,0,0};
	int matrx_mlt2_test_input_1_m = 3;
	int matrx_mlt2_test_input_1_n = 3;
	int matrx_mlt2_test_input_b_1_m = 3;
	int matrx_mlt2_test_input_b_1_n = 1;
	test_matrx_mlt2(matrx_mlt2_out_1,matrx_mlt2_test_input_1,matrx_mlt2_test_input_1_m,matrx_mlt2_test_input_1_n,matrx_mlt2_test_input_b_1,matrx_mlt2_test_input_b_1_m,matrx_mlt2_test_input_b_1_n,matrx_mlt2_expect_1,epsilon);

	//matrx_sub3

	double sub3_test_input_1[]={1,2,3,4};
	double sub3_test_input_d_1=5;
	double sub3_expect_1[]={-4,-3,-2,-1};
	int sub3_length_on_side_1= sqrt(sizeof(sub3_test_input_1) / sizeof(double));
	test_matrx_sub3(sub3_test_input_1,sub3_expect_1,sub3_test_input_d_1,sub3_length_on_side_1,sub3_length_on_side_1,epsilon);

	//matrx_det

	double det_test_input_1[]={1.1776,0.7765,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0.7765,1.1776,0,0.1629,0.3000,0.1629,0.4960,0.0670,0.0670,0.0670,0.4960,0.0670,0,0,1.1776,0,0,0,0,0,0,0,0,0,0.1629,0.1629,0,1.1776,0.1629,0.3000,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.3000,0.3000,0,0.1629,1.1776,0.1629,0.3000,0.0670,0.0670,0.0670,0.3000,0.0670,0.1629,0.1629,0,0.3000,0.1629,1.1776,0.1629,0.0670,0.0670,0.0670,0.1629,0.0670,0.4960,0.4960,0,0.1629,0.3000,0.1629,1.1776,0.0670,0.0670,0.0670,0.7765,0.0670,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,1.1776,0.3000,0.3000,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.3000,1.1776,0.4960,0.0670,0.1629,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.3000,0.4960,1.1776,0.0670,0.1629,0.4960,0.4960,0,0.1629,0.3000,0.1629,0.7765,0.0670,0.0670,0.0670,1.1776,0.0670,0.0670,0.0670,0,0.0670,0.0670,0.0670,0.0670,0.1629,0.1629,0.1629,0.0670,1.1776};
	double expect_1=0.9351;
	int length_on_side_1 = sqrt(sizeof(det_test_input_1) / sizeof(double));
	test_matrx_det(det_test_input_1,expect_1,length_on_side_1,big_epsilon); 

	double det_test_input_2[]={1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1};
	double expect_2=1;
	int length_on_side_2 = sqrt(sizeof(det_test_input_2) / sizeof(double));
	test_matrx_det(det_test_input_2,expect_2,length_on_side_2,big_epsilon);

	double det_test_input_3[]={NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	double expect_3=NAN;
	int length_on_side_3 = sqrt(sizeof(det_test_input_3) / sizeof(double));
	test_matrx_det(det_test_input_3,expect_3,length_on_side_3,big_epsilon);

	double det_test_input_4[]={1.3482,1.0327,0,0.3243,0.5263,0.3243,0.7607,0.1501,0.1501,0.1501,0.7607,0.1501,1.0327,1.3482,0,0.3243,0.5263,0.3243,0.7607,0.1501,0.1501,0.1501,0.7607,0.1501,0,0,1.3482,0,0,0,0,0,0,0,0,0,0.3243,0.3243,0,1.3482,0.3243,0.5263,0.3243,0.1501,0.1501,0.1501,0.3243,0.1501,0.5263,0.5263,0,0.3243,1.3482,0.3243,0.5263,0.1501,0.1501,0.1501,0.5263,0.1501,0.3243,0.3243,0,0.5263,0.3243,1.3482,0.3243,0.1501,0.1501,0.1501,0.3243,0.1501,0.7607,0.7607,0,0.3243,0.5263,0.3243,1.3482,0.1501,0.1501,0.1501,1.0327,0.1501,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,1.3482,0.5263,0.5263,0.1501,0.3243,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,0.5263,1.3482,0.7607,0.1501,0.3243,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,0.5263,0.7607,1.3482,0.1501,0.3243,0.7607,0.7607,0,0.3243,0.5263,0.3243,1.0327,0.1501,0.1501,0.1501,1.3482,0.1501,0.1501,0.1501,0,0.1501,0.1501,0.1501,0.1501,0.3243,0.3243,0.3243,0.1501,1.3482};
	double expect_4=0.9884;
	int length_on_side_4 = sqrt(sizeof(det_test_input_4) / sizeof(double));
	test_matrx_det(det_test_input_4,expect_4,length_on_side_4,big_epsilon);

	double det_test_input_5[]={1.0521,0.5069,0,0.0464,0.1100,0.0464,0.2403,0.0152,0.0152,0.0152,0.2403,0.0152,0.5069,1.0521,0,0.0464,0.1100,0.0464,0.2403,0.0152,0.0152,0.0152,0.2403,0.0152,0,0,1.0521,0,0,0,0,0,0,0,0,0,0.0464,0.0464,0,1.0521,0.0464,0.1100,0.0464,0.0152,0.0152,0.0152,0.0464,0.0152,0.1100,0.1100,0,0.0464,1.0521,0.0464,0.1100,0.0152,0.0152,0.0152,0.1100,0.0152,0.0464,0.0464,0,0.1100,0.0464,1.0521,0.0464,0.0152,0.0152,0.0152,0.0464,0.0152,0.2403,0.2403,0,0.0464,0.1100,0.0464,1.0521,0.0152,0.0152,0.0152,0.5069,0.0152,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,1.0521,0.1100,0.1100,0.0152,0.0464,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,0.1100,1.0521,0.2403,0.0152,0.0464,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,0.1100,0.2403,1.0521,0.0152,0.0464,0.2403,0.2403,0,0.0464,0.1100,0.0464,0.5069,0.0152,0.0152,0.0152,1.0521,0.0152,0.0152,0.0152,0,0.0152,0.0152,0.0152,0.0152,0.0464,0.0464,0.0464,0.0152,1.0521};
	double expect_5=0.8688;
	int length_on_side_5 = sqrt(sizeof(det_test_input_5) / sizeof(double));
	test_matrx_det(det_test_input_5,expect_5,length_on_side_5,big_epsilon);

	double det_test_input_6[]={1.1941,0.1856,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1856,1.1941,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,1.1941,0.1856,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,0.1856,1.1941,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.1941,0.6187,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.6187,1.1941,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,1.1941,0.6187,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,0.6187,1.1941,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,1.1941,0.2916,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.2916,1.1941,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,1.1941,0.8661,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.8661,1.1941,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.6187,0.6187,1.1941,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,1.1941,0.6187,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,1.1941,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,0.6187,1.1941,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0.0454,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0.0454,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.0454,0.0454,1.1941};
	double expect_6=0.9226;
	int length_on_side_6 = sqrt(sizeof(det_test_input_6) / sizeof(double));
	test_matrx_det(det_test_input_6,expect_6,length_on_side_6,big_epsilon);

	double det_test_input_7[]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	double expect_7=1;
	int length_on_side_7 = sqrt(sizeof(det_test_input_7) / sizeof(double));
	test_matrx_det(det_test_input_7,expect_7,length_on_side_7,big_epsilon);

	double det_test_input_8[]={NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	double expect_8=NAN;
	int length_on_side_8 = sqrt(sizeof(det_test_input_8) / sizeof(double));
	test_matrx_det(det_test_input_8,expect_8,length_on_side_8,big_epsilon);

	double det_test_input_9[]={1.3908,0.3771,0.2364,0.2364,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.3771,1.3908,0.2364,0.2364,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2364,0.2364,1.3908,0.3771,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2364,0.2364,0.3771,1.3908,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,1.3908,0.2364,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,0.2364,1.3908,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.3908,0.9131,0.7132,0.7132,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.9131,1.3908,0.7132,0.7132,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.7132,0.7132,1.3908,0.9131,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.7132,0.7132,0.9131,1.3908,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.5353,0.5353,0.5353,0.5353,1.3908,0.9131,0.7132,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.5353,0.5353,0.5353,0.5353,0.9131,1.3908,0.7132,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.5353,0.5353,0.5353,0.5353,0.7132,0.7132,1.3908,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,1.3908,0.5353,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.5353,1.3908,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,1.3908,1.1380,0.9131,0.7132,0.7132,0.7132,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,1.1380,1.3908,0.9131,0.7132,0.7132,0.7132,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.9131,0.9131,1.3908,0.7132,0.7132,0.7132,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.7132,0.7132,0.7132,1.3908,0.9131,0.9131,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.7132,0.7132,0.7132,0.9131,1.3908,0.9131,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.7132,0.7132,0.7132,0.9131,0.9131,1.3908,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.5353,0.5353,0.5353,0.5353,0.5353,0.5353,1.3908,0.9131,0.7132,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.5353,0.5353,0.5353,0.5353,0.5353,0.5353,0.9131,1.3908,0.7132,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.5353,0.5353,0.5353,0.5353,0.5353,0.5353,0.7132,0.7132,1.3908,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,1.3908,0.2364,0.1113,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.2364,1.3908,0.1113,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.1113,0.1113,1.3908};
	double expect_9=0.9975;
	int length_on_side_9 = sqrt(sizeof(det_test_input_9) / sizeof(double));
	test_matrx_det(det_test_input_9,expect_9,length_on_side_9,big_epsilon);

	double det_test_input_10[]={1.0550,0.0517,0.0244,0.0244,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0517,1.0550,0.0244,0.0244,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0244,0.0244,1.0550,0.0517,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0244,0.0244,0.0517,1.0550,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,1.0550,0.0244,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,0.0244,1.0550,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0550,0.3331,0.1843,0.1843,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.3331,1.0550,0.1843,0.1843,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.1843,0.1843,1.0550,0.3331,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.1843,0.1843,0.3331,1.0550,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0998,0.0998,0.0998,0.0998,1.0550,0.3331,0.1843,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0998,0.0998,0.0998,0.0998,0.3331,1.0550,0.1843,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0998,0.0998,0.0998,0.0998,0.1843,0.1843,1.0550,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,1.0550,0.0998,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0998,1.0550,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,1.0550,0.5947,0.3331,0.1843,0.1843,0.1843,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.5947,1.0550,0.3331,0.1843,0.1843,0.1843,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.3331,0.3331,1.0550,0.1843,0.1843,0.1843,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.1843,0.1843,0.1843,1.0550,0.3331,0.3331,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.1843,0.1843,0.1843,0.3331,1.0550,0.3331,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.1843,0.1843,0.1843,0.3331,0.3331,1.0550,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0998,0.0998,0.0998,0.0998,0.0998,0.0998,1.0550,0.3331,0.1843,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0998,0.0998,0.0998,0.0998,0.0998,0.0998,0.3331,1.0550,0.1843,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0998,0.0998,0.0998,0.0998,0.0998,0.0998,0.1843,0.1843,1.0550,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,1.0550,0.0244,0.0088,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0244,1.0550,0.0088,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0088,0.0088,1.0550};
	double expect_10=0.8330;
	int length_on_side_10 = sqrt(sizeof(det_test_input_10) / sizeof(double));
	test_matrx_det(det_test_input_10,expect_10,length_on_side_10,big_epsilon);

	double det_test_input_11[]={1,2,3,4};
	double expect_11=-2;
	int length_on_side_11 = sqrt(sizeof(det_test_input_11) / sizeof(double));
	test_matrx_det(det_test_input_11,expect_11,length_on_side_11,big_epsilon);

	double det_test_input_12[]={3,0,6,1};
	double expect_12=3;
	int length_on_side_12 = sqrt(sizeof(det_test_input_12) / sizeof(double));
	test_matrx_det(det_test_input_12,expect_12,length_on_side_12,big_epsilon);

	double det_test_input_13[]={3,0,0,1};
	double expect_13=3;
	int length_on_side_13 = sqrt(sizeof(det_test_input_13) / sizeof(double));
	test_matrx_det(det_test_input_13,expect_13,length_on_side_13,big_epsilon);

	double det_test_input_14[]={3,0,4,0,1,0,4,0,0};
	double expect_14=-16;
	int length_on_side_14 = sqrt(sizeof(det_test_input_14) / sizeof(double));
	test_matrx_det(det_test_input_14,expect_14,length_on_side_14,big_epsilon);

	//test funct steps with several values

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

	double expect_vp_5_5[]= {1.1941,0.1856,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1856,1.1941,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,1.1941,0.1856,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,0.1856,1.1941,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.1941,0.6187,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.6187,1.1941,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,1.1941,0.6187,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,0.6187,1.1941,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,1.1941,0.2916,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.2916,1.1941,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,1.1941,0.8661,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.8661,1.1941,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.6187,0.6187,1.1941,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,1.1941,0.6187,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,1.1941,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,0.6187,1.1941,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0.0454,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0.0454,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.0454,0.0454,1.1941};
	int expect_vp_5_5_size = sizeof(expect_vp_5_5) / sizeof(double);
	test_vp_calculation(0.5,0.5,expect_vp_5_5,expect_vp_5_5_size,epsilon);

	double expect_vp_0_5[]= {1.1941,0.1856,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1856,1.1941,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,1.1941,0.1856,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,0.1856,1.1941,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.1941,0.6187,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.6187,1.1941,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,1.1941,0.6187,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,0.6187,1.1941,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,1.1941,0.2916,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.2916,1.1941,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,1.1941,0.8661,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.8661,1.1941,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.6187,0.6187,1.1941,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,1.1941,0.6187,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,1.1941,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,0.6187,1.1941,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0.0454,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0.0454,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.0454,0.0454,1.1941};
	int expect_vp_0_5_size = sizeof(expect_vp_0_5) / sizeof(double);
	test_vp_calculation(0.0,0.5,expect_vp_0_5,expect_vp_0_5_size,epsilon);

	double expect_vp_5_0[]= {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	int expect_vp_5_0_size = sizeof(expect_vp_5_0) / sizeof(double);
	test_vp_calculation(0.5,0.0,expect_vp_5_0,expect_vp_5_0_size,epsilon);

	double expect_vp_1_5[]= {1.1941,0.1856,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1856,1.1941,0.1057,0.1057,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,1.1941,0.1856,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1057,0.1057,0.1856,1.1941,0.0454,0.0454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.1941,0.6187,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.6187,1.1941,0.4322,0.4322,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,1.1941,0.6187,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.4322,0.4322,0.6187,1.1941,0.2916,0.2916,0.2916,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.1856,0.1856,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,1.1941,0.2916,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.2916,1.1941,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.0454,0.0454,0.1057,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,1.1941,0.8661,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.8661,1.1941,0.6187,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.6187,0.6187,1.1941,0.4322,0.4322,0.4322,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,1.1941,0.6187,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,1.1941,0.6187,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.4322,0.4322,0.4322,0.6187,0.6187,1.1941,0.2916,0.2916,0.2916,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,1.1941,0.6187,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.6187,1.1941,0.4322,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.2916,0.2916,0.2916,0.2916,0.2916,0.2916,0.4322,0.4322,1.1941,0.0454,0.0454,0.1856,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,1.1941,0.1057,0.0454,0,0,0,0,0,0,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.0454,0.1057,1.1941,0.0454,0,0,0,0,0,0,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1057,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.1856,0.0454,0.0454,1.1941};
	int expect_vp_1_5_size = sizeof(expect_vp_1_5) / sizeof(double);
	test_vp_calculation(1.0,0.5,expect_vp_1_5,expect_vp_1_5_size,epsilon);

	double expect_vp_5_1[]= {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	int expect_vp_5_1_size = sizeof(expect_vp_5_1) / sizeof(double);
	test_vp_calculation(0.5,1.0,expect_vp_5_1,expect_vp_5_1_size,epsilon);

	double expect_vp_0_0[]= {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
	int expect_vp_0_0_size = sizeof(expect_vp_0_0) / sizeof(double);
	test_vp_calculation(0.0,0.0,expect_vp_0_0,expect_vp_0_0_size,epsilon);

	double expect_vp_1_1[]= {NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	int expect_vp_1_1_size = sizeof(expect_vp_1_1) / sizeof(double);
	test_vp_calculation(1.0,1.0,expect_vp_1_1,expect_vp_1_1_size,epsilon);

	double expect_vp_75_75[]= {1.3908,0.3771,0.2364,0.2364,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.3771,1.3908,0.2364,0.2364,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2364,0.2364,1.3908,0.3771,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2364,0.2364,0.3771,1.3908,0.1113,0.1113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,1.3908,0.2364,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,0.2364,1.3908,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.3908,0.9131,0.7132,0.7132,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.9131,1.3908,0.7132,0.7132,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.7132,0.7132,1.3908,0.9131,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.7132,0.7132,0.9131,1.3908,0.5353,0.5353,0.5353,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.5353,0.5353,0.5353,0.5353,1.3908,0.9131,0.7132,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.5353,0.5353,0.5353,0.5353,0.9131,1.3908,0.7132,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.5353,0.5353,0.5353,0.5353,0.7132,0.7132,1.3908,0.3771,0.3771,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,1.3908,0.5353,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.5353,1.3908,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.1113,0.1113,0.2364,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,1.3908,1.1380,0.9131,0.7132,0.7132,0.7132,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,1.1380,1.3908,0.9131,0.7132,0.7132,0.7132,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.9131,0.9131,1.3908,0.7132,0.7132,0.7132,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.7132,0.7132,0.7132,1.3908,0.9131,0.9131,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.7132,0.7132,0.7132,0.9131,1.3908,0.9131,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.7132,0.7132,0.7132,0.9131,0.9131,1.3908,0.5353,0.5353,0.5353,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.5353,0.5353,0.5353,0.5353,0.5353,0.5353,1.3908,0.9131,0.7132,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.5353,0.5353,0.5353,0.5353,0.5353,0.5353,0.9131,1.3908,0.7132,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.5353,0.5353,0.5353,0.5353,0.5353,0.5353,0.7132,0.7132,1.3908,0.1113,0.1113,0.3771,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,1.3908,0.2364,0.1113,0,0,0,0,0,0,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.1113,0.2364,1.3908,0.1113,0,0,0,0,0,0,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.2364,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.3771,0.1113,0.1113,1.3908};
	int expect_vp_75_75_size = sizeof(expect_vp_75_75) / sizeof(double);
	test_vp_calculation(0.75,0.75,expect_vp_75_75,expect_vp_75_75_size,epsilon);

	double expect_vp_25_25[]= {1.0550,0.0517,0.0244,0.0244,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0517,1.0550,0.0244,0.0244,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0244,0.0244,1.0550,0.0517,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0244,0.0244,0.0517,1.0550,0.0088,0.0088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,1.0550,0.0244,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,0.0244,1.0550,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0550,0.3331,0.1843,0.1843,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.3331,1.0550,0.1843,0.1843,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.1843,0.1843,1.0550,0.3331,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.1843,0.1843,0.3331,1.0550,0.0998,0.0998,0.0998,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0998,0.0998,0.0998,0.0998,1.0550,0.3331,0.1843,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0998,0.0998,0.0998,0.0998,0.3331,1.0550,0.1843,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0998,0.0998,0.0998,0.0998,0.1843,0.1843,1.0550,0.0517,0.0517,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,1.0550,0.0998,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0998,1.0550,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0088,0.0088,0.0244,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,1.0550,0.5947,0.3331,0.1843,0.1843,0.1843,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.5947,1.0550,0.3331,0.1843,0.1843,0.1843,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.3331,0.3331,1.0550,0.1843,0.1843,0.1843,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.1843,0.1843,0.1843,1.0550,0.3331,0.3331,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.1843,0.1843,0.1843,0.3331,1.0550,0.3331,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.1843,0.1843,0.1843,0.3331,0.3331,1.0550,0.0998,0.0998,0.0998,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0998,0.0998,0.0998,0.0998,0.0998,0.0998,1.0550,0.3331,0.1843,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0998,0.0998,0.0998,0.0998,0.0998,0.0998,0.3331,1.0550,0.1843,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0998,0.0998,0.0998,0.0998,0.0998,0.0998,0.1843,0.1843,1.0550,0.0088,0.0088,0.0517,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,1.0550,0.0244,0.0088,0,0,0,0,0,0,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0088,0.0244,1.0550,0.0088,0,0,0,0,0,0,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0244,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0517,0.0088,0.0088,1.0550};
	int expect_vp_25_25_size = sizeof(expect_vp_25_25) / sizeof(double);
	test_vp_calculation(0.25,0.25,expect_vp_25_25,expect_vp_25_25_size,epsilon);

	/*Vh=Vh./det(Vh)^(1/p);-----------------------------------*/

	double expect_vh_rdv_5_5[]={1.1842,0.7808,0,0.1638,0.3016,0.1638,0.4988,0.0674,0.0674,0.0674,0.4988,0.0674,0.7808,1.1842,0,0.1638,0.3016,0.1638,0.4988,0.0674,0.0674,0.0674,0.4988,0.0674,0,0,1.1842,0,0,0,0,0,0,0,0,0,0.1638,0.1638,0,1.1842,0.1638,0.3016,0.1638,0.0674,0.0674,0.0674,0.1638,0.0674,0.3016,0.3016,0,0.1638,1.1842,0.1638,0.3016,0.0674,0.0674,0.0674,0.3016,0.0674,0.1638,0.1638,0,0.3016,0.1638,1.1842,0.1638,0.0674,0.0674,0.0674,0.1638,0.0674,0.4988,0.4988,0,0.1638,0.3016,0.1638,1.1842,0.0674,0.0674,0.0674,0.7808,0.0674,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,1.1842,0.3016,0.3016,0.0674,0.1638,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.3016,1.1842,0.4988,0.0674,0.1638,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.3016,0.4988,1.1842,0.0674,0.1638,0.4988,0.4988,0,0.1638,0.3016,0.1638,0.7808,0.0674,0.0674,0.0674,1.1842,0.0674,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.1638,0.1638,0.1638,0.0674,1.1842};
	int expect_vh_rdv_5_5_size = sizeof(expect_vh_rdv_5_5) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_5_5,expect_vh_rdv_5_5,expect_vh_rdv_5_5_size,big_epsilon);

	double expect_vh_rdv_0_5[]={1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1};
	int expect_vh_rdv_0_5_size = sizeof(expect_vh_rdv_0_5) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_0_5,expect_vh_rdv_0_5,expect_vh_rdv_0_5_size,big_epsilon);

	double expect_vh_rdv_5_0[]={1.1842,0.7808,0,0.1638,0.3016,0.1638,0.4988,0.0674,0.0674,0.0674,0.4988,0.0674,0.7808,1.1842,0,0.1638,0.3016,0.1638,0.4988,0.0674,0.0674,0.0674,0.4988,0.0674,0,0,1.1842,0,0,0,0,0,0,0,0,0,0.1638,0.1638,0,1.1842,0.1638,0.3016,0.1638,0.0674,0.0674,0.0674,0.1638,0.0674,0.3016,0.3016,0,0.1638,1.1842,0.1638,0.3016,0.0674,0.0674,0.0674,0.3016,0.0674,0.1638,0.1638,0,0.3016,0.1638,1.1842,0.1638,0.0674,0.0674,0.0674,0.1638,0.0674,0.4988,0.4988,0,0.1638,0.3016,0.1638,1.1842,0.0674,0.0674,0.0674,0.7808,0.0674,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,1.1842,0.3016,0.3016,0.0674,0.1638,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.3016,1.1842,0.4988,0.0674,0.1638,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.3016,0.4988,1.1842,0.0674,0.1638,0.4988,0.4988,0,0.1638,0.3016,0.1638,0.7808,0.0674,0.0674,0.0674,1.1842,0.0674,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.1638,0.1638,0.1638,0.0674,1.1842};
	int expect_vh_rdv_5_0_size = sizeof(expect_vh_rdv_5_0) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_5_0,expect_vh_rdv_5_0,expect_vh_rdv_5_0_size,big_epsilon);

	double expect_vh_rdv_1_5[]={NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	int expect_vh_rdv_1_5_size = sizeof(expect_vh_rdv_1_5) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_1_5,expect_vh_rdv_1_5,expect_vh_rdv_1_5_size,big_epsilon);

	double expect_vh_rdv_5_1[]={1.1842,0.7808,0,0.1638,0.3016,0.1638,0.4988,0.0674,0.0674,0.0674,0.4988,0.0674,0.7808,1.1842,0,0.1638,0.3016,0.1638,0.4988,0.0674,0.0674,0.0674,0.4988,0.0674,0,0,1.1842,0,0,0,0,0,0,0,0,0,0.1638,0.1638,0,1.1842,0.1638,0.3016,0.1638,0.0674,0.0674,0.0674,0.1638,0.0674,0.3016,0.3016,0,0.1638,1.1842,0.1638,0.3016,0.0674,0.0674,0.0674,0.3016,0.0674,0.1638,0.1638,0,0.3016,0.1638,1.1842,0.1638,0.0674,0.0674,0.0674,0.1638,0.0674,0.4988,0.4988,0,0.1638,0.3016,0.1638,1.1842,0.0674,0.0674,0.0674,0.7808,0.0674,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,1.1842,0.3016,0.3016,0.0674,0.1638,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.3016,1.1842,0.4988,0.0674,0.1638,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.3016,0.4988,1.1842,0.0674,0.1638,0.4988,0.4988,0,0.1638,0.3016,0.1638,0.7808,0.0674,0.0674,0.0674,1.1842,0.0674,0.0674,0.0674,0,0.0674,0.0674,0.0674,0.0674,0.1638,0.1638,0.1638,0.0674,1.1842};
	int expect_vh_rdv_5_1_size = sizeof(expect_vh_rdv_5_1) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_5_1,expect_vh_rdv_5_1,expect_vh_rdv_5_1_size,big_epsilon);

	double expect_vh_rdv_0_0[]={1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1};
	int expect_vh_rdv_0_0_size = sizeof(expect_vh_rdv_0_0) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_0_0,expect_vh_rdv_0_0,expect_vh_rdv_0_0_size,big_epsilon);

	double expect_vh_rdv_1_1[]={NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
	int expect_vh_rdv_1_1_size = sizeof(expect_vh_rdv_1_1) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_1_1,expect_vh_rdv_1_1,expect_vh_rdv_1_1_size,big_epsilon);

	double expect_vh_rdv_75_75[]={1.3495,1.0337,0,0.3246,0.5268,0.3246,0.7615,0.1503,0.1503,0.1503,0.7615,0.1503,1.0337,1.3495,0,0.3246,0.5268,0.3246,0.7615,0.1503,0.1503,0.1503,0.7615,0.1503,0,0,1.3495,0,0,0,0,0,0,0,0,0,0.3246,0.3246,0,1.3495,0.3246,0.5268,0.3246,0.1503,0.1503,0.1503,0.3246,0.1503,0.5268,0.5268,0,0.3246,1.3495,0.3246,0.5268,0.1503,0.1503,0.1503,0.5268,0.1503,0.3246,0.3246,0,0.5268,0.3246,1.3495,0.3246,0.1503,0.1503,0.1503,0.3246,0.1503,0.7615,0.7615,0,0.3246,0.5268,0.3246,1.3495,0.1503,0.1503,0.1503,1.0337,0.1503,0.1503,0.1503,0,0.1503,0.1503,0.1503,0.1503,1.3495,0.5268,0.5268,0.1503,0.3246,0.1503,0.1503,0,0.1503,0.1503,0.1503,0.1503,0.5268,1.3495,0.7615,0.1503,0.3246,0.1503,0.1503,0,0.1503,0.1503,0.1503,0.1503,0.5268,0.7615,1.3495,0.1503,0.3246,0.7615,0.7615,0,0.3246,0.5268,0.3246,1.0337,0.1503,0.1503,0.1503,1.3495,0.1503,0.1503,0.1503,0,0.1503,0.1503,0.1503,0.1503,0.3246,0.3246,0.3246,0.1503,1.3495};
	int expect_vh_rdv_75_75_size = sizeof(expect_vh_rdv_75_75) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_75_75,expect_vh_rdv_75_75,expect_vh_rdv_75_75_size,big_epsilon);

	double expect_vh_rdv_25_25[]={1.0645,0.5128,0,0.0469,0.1113,0.0469,0.2432,0.0154,0.0154,0.0154,0.2432,0.0154,0.5128,1.0645,0,0.0469,0.1113,0.0469,0.2432,0.0154,0.0154,0.0154,0.2432,0.0154,0,0,1.0645,0,0,0,0,0,0,0,0,0,0.0469,0.0469,0,1.0645,0.0469,0.1113,0.0469,0.0154,0.0154,0.0154,0.0469,0.0154,0.1113,0.1113,0,0.0469,1.0645,0.0469,0.1113,0.0154,0.0154,0.0154,0.1113,0.0154,0.0469,0.0469,0,0.1113,0.0469,1.0645,0.0469,0.0154,0.0154,0.0154,0.0469,0.0154,0.2432,0.2432,0,0.0469,0.1113,0.0469,1.0645,0.0154,0.0154,0.0154,0.5128,0.0154,0.0154,0.0154,0,0.0154,0.0154,0.0154,0.0154,1.0645,0.1113,0.1113,0.0154,0.0469,0.0154,0.0154,0,0.0154,0.0154,0.0154,0.0154,0.1113,1.0645,0.2432,0.0154,0.0469,0.0154,0.0154,0,0.0154,0.0154,0.0154,0.0154,0.1113,0.2432,1.0645,0.0154,0.0469,0.2432,0.2432,0,0.0469,0.1113,0.0469,0.5128,0.0154,0.0154,0.0154,1.0645,0.0154,0.0154,0.0154,0,0.0154,0.0154,0.0154,0.0154,0.0469,0.0469,0.0469,0.0154,1.0645};
	int expect_vh_rdv_25_25_size = sizeof(expect_vh_rdv_25_25) / sizeof(double);
	test_vh_rdv_calculation(expect_vh_25_25,expect_vh_rdv_25_25,expect_vh_rdv_25_25_size,big_epsilon);

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

	printf("est=  %f  %f\n",*(XMIN),*(XMIN+1));
	printf("MSE=  %f\n",YNEWLO);
	if(all_passed == 1) {
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
	printf(" -> test_funct");
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

	if(vh_size == correct_size) {
		for(int i=0; i <correct_size; i++) {
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
		printf(" -> %i correct\n",correct_size);
	}
}

void test_vp_calculation(double d1, double d2, double correct_return[], int correct_size, double epsilon) {

	double A_QQ[q*q];
	matrx_mlt(A_QQ,2.00,initVp,q,q);
	array_pow(A_QQ,d2,A_QQ,q,q);
	matrx_sub(A_QQ,1.00,A_QQ,q,q);

	double B_QQ[q*q];
	array_pow(B_QQ,d2,tau2,q,q);

	array_mlt(B_QQ,B_QQ,q,q,A_QQ);
	array_rdv(B_QQ,B_QQ,q,q,1.00 - d2*d2);

	double vp_size = sizeof(B_QQ) / sizeof(double);

	if(vp_size == correct_size) {
		for(int i=0; i <correct_size; i++) {
			double actual = B_QQ[i];
			double expect = correct_return[i];
			double error = actual - expect;
			if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
				printf("* PASSED");
			} else {
				all_passed = 0;
				printf("X FAILED");
			}
			printf(" -> test_vp_calculation at %i",i);
			printf(" -> d1 = %f",d1);
			printf(" -> d2 = %f",d2);
			printf(" -> %f returned",actual);
			printf(" -> %f correct\n",expect);
		}
	} else {
		all_passed = 0;
		printf("X FAILED");
		printf(" -> test_vp_calculation() -- size");
		printf(" -> d1 = %f",d1);
		printf(" -> d2 = %f",d2);
		printf(" -> %f returned",vp_size);
		printf(" -> %i correct\n",correct_size);
	}
}

void test_matrx_det(double A[], double expect, int size, double epsilon){

	double actual = matrx_det(A,size);

	double error = actual - expect;

	if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
		printf("* PASSED");
	} else {
		all_passed = 0;
		printf("X FAILED");
	}
	printf(" -> test_matrx_det");
	printf(" -> %f returned",actual);
	printf(" -> %f correct\n",expect);
}

void test_vh_rdv_calculation(double io_vh[], double correct_return[], int correct_size, double epsilon){

	double detrm = matrx_det(io_vh,p);
	double odp = (1.00/(double) p);

	printf("test_vh_rdv matrx_det is %f ",detrm);

	array_rdv(io_vh,io_vh,p,p,pow(detrm,odp));

	for(int i=0; i <correct_size; i++) {
		double actual = io_vh[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_vh_rdv_calculation at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void output(double *matrix,int m,int n) {

	int i;
	int j;
	for(i=0; i<m; i++) {
		for(j=0; j<n; j++) {
			printf("%f ",*(matrix+(i*m+j)));
		}
		printf("\n");
	}
}

void test_matrx_mlt(double io[], double correct_return[], double d, double m, double n, double epsilon){

	matrx_mlt(io,d,io,m,n);

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_matrx_mlt at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_array_pow(double io[],double correct_return[], double d,double m,double n,double epsilon){

	array_pow(io,d,io,m,n);

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_array_pow at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_matrx_sub(double io[],double correct_return[], double d,double m,double n,double epsilon){

	matrx_sub(io,d,io,m,n);

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_matrx_sub at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_array_amlt(double io[],double correct_return[],double b[],double m,double n,double epsilon){

	array_mlt(io,io,m,n,b);

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_array_amlt at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_array_rdv(double io[],double correct_return[],double d,double m,double n,double epsilon){

	array_rdv(io,io,m,n,d);

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_array_rdv at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_kron(double io[],double correct_return[],double b[],double ma,double na,double mb,double nb,double epsilon){

	kron(io,io,ma,na,b,mb,nb);

	int m = ma * mb;
	int n = na * nb;

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_kron at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_tran(double io[],double correct_return[],double m,double n,double epsilon){

	tran(io,io,m,n);

	for(int i=0; i <(m * n); i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_tran at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_matrx_inv_(double io[],double correct_return[],int n,double epsilon){

	matrx_inv(io,io,n);

	int m = n;

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_matrx_inv at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_ones(double io[],double correct_return[],int m,int n,double epsilon){

	ones(io,m,n);

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_ones at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_matrx_mlt2(double io[],double io1[],int m,int n,double io2[],int ma,int na,double correct_return[],double epsilon){

	matrx_mlt2(io,io1,m,n,io2,ma,na);

	for(int i=0; i <m * na; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_matrx_mlt2 at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}

void test_matrx_sub3(double io[],double correct_return[], double d,double m,double n,double epsilon){

	matrx_sub3(io,io,m,n,d);

	for(int i=0; i <m * n; i++) {
		double actual = io[i];
		double expect = correct_return[i];
		double error = actual - expect;
		if( isnan(actual) && isnan(expect) || (error * error) < epsilon ) {
			printf("* PASSED");
		} else {
			all_passed = 0;
			printf("X FAILED");
		}
		printf(" -> test_matrx_sub3 at %i",i);
		printf(" -> %f returned",actual);
		printf(" -> %f correct\n",expect);
	}
}
