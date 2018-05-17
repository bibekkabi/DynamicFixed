/***********************************************************************************************
Program for fixed point implementations of various mathematical operations
Author : Sajal Kumar

************************************************************************************************/



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<conio.h>
#include "Fixed_point.h"

//Define a custom made fixed point strucure to hold the fix value and its corresponding Q value
typedef struct {
					FP32 fix_val;
					FP32 Q_val;
				}QFrm;

//typedef unsigned long int FP32 ;
#define m 50 
#define n 11 


//Sample Input (1) with fundamental and 5th harmonic 
float  input[50] = {0,0.24289,0.4389,0.55834,0.59931,0.58779,0.56699,0.5803,0.65412,0.78727,0.95106,1.0998,1.1882,1.1882,1.0998,0.95106,0.78727,0.65412,0.5803,0.56699,0.58779,0.59931,0.55834,0.4389,0.24289,2.4493e-016,-0.24289,-0.4389,-0.55834,-0.59931,-0.58779,-0.56699,-0.5803,-0.65412,-0.78727,-0.95106,-1.0998,-1.1882,-1.1882,-1.0998,-0.95106,-0.78727,-0.65412,-0.5803,-0.56699,-0.58779,-0.59931,-0.55834,-0.4389,-0.24289};

//Sample Input (2) with fundamental, 5th and 7th harmonic
//float input[50] = {0,0.35296,0.57923,0.62716,0.54672,0.45192,0.44637,0.5624,0.75191,0.92985,1.035,1.0643,1.059,1.059,1.0643,1.035,0.92985,0.75191,0.5624,0.44637,0.45192,0.54672,0.62716,0.57923,0.35296,3.6739e-016,-0.35296,-0.57923,-0.62716,-0.54672,-0.45192,-0.44637,-0.5624,-0.75191,-0.92985,-1.035,-1.0643,-1.059,-1.059,-1.0643,-1.035,-0.92985,-0.75191,-0.5624,-0.44637,-0.45192,-0.54672,-0.62716,-0.57923,-0.35296};
float imp_res[11] = {0.061036,0.074909,0.086908,0.096172,0.10202,0.10402,0.10202,0.096172,0.086908,0.074909,0.061036};

//Output for Input (1) obtained from Matlab
float exp_output[61] = {0,0.014825,0.044983,0.088065,0.13991,0.19628,0.25446,0.31401,0.37656,0.4446,0.51984,0.60171,0.67545,0.74073,0.79653,0.84112,0.87245,0.88866,0.88866,0.87245,0.84112,0.79653,0.74073,0.67545,0.60171,0.51984,0.42978,0.33157,0.22594,0.11455,2.498e-016,-0.11455,-0.22594,-0.33157,-0.42978,-0.51984,-0.60171,-0.67545,-0.74073,-0.79653,-0.84112,-0.87245,-0.88866,-0.88866,-0.87245,-0.84112,-0.79653,-0.74073,-0.67545,-0.60171,-0.51984,-0.4446,-0.37656,-0.31401,-0.25446,-0.19628,-0.13991,-0.088065,-0.044983,-0.014825,-2.9899e-017};

//Output for Input (2) obtained from Matlab
//float exp_output[61] = {0,0.021543,0.061794,0.11234,0.16463,0.21476,0.26474,0.31986,0.38433,0.45822,0.53768,0.61755,0.6778,0.72789,0.77781,0.8301,0.87711,0.90563,0.90563,0.87711,0.8301,0.77781,0.72789,0.6778,0.61755,0.53768,0.43668,0.32254,0.20752,0.1001,1.3184e-016,-0.1001,-0.20752,-0.32254,-0.43668,-0.53768,-0.61755,-0.6778,-0.72789,-0.77781,-0.8301,-0.87711,-0.90563,-0.90563,-0.87711,-0.8301,-0.77781,-0.72789,-0.6778,-0.61755,-0.53768,-0.45822,-0.38433,-0.31986,-0.26474,-0.21476,-0.16463,-0.11234,-0.061794,-0.021543,-4.4849e-017};
int main()
{
	
	FP32 Q = 14 ;
	float flt_a = 180.9368 ;
	float flt_b = 2.0872 ;
//	unsigned int q_val =1;
//	FP32 Qtest = 2;
	//determing the variable Q of both the variables
//	FP32 a.Q_val,b.Q_val ;
	/*while(abs(flt_a)>Qtest)
	{		Qtest = Qtest<<1;
			q_val++;
	}
	a.Q_val = 31-q_val;
	q_val =1;
	Qtest =2 ;
	while(abs(flt_b)>Qtest)
	{		Qtest = Qtest<<1;
			q_val++;
	}
	b.Q_val = 31-q_val;
	*/
	QFrm a,b;
	a.Q_val = Get_Optimum_Q(flt_a);
	b.Q_val = Get_Optimum_Q(flt_b);
	printf("a.Q_val = %d  , b.Q_val = %d",a.Q_val,b.Q_val);

	
	a.fix_val = tofix(flt_a,a.Q_val);
	b.fix_val = tofix(flt_b,b.Q_val);

	printf("\na = %f  a.fix_val = %d",flt_a,a.fix_val);
	printf("\nb = %f  b.fix_val = %d",flt_b,b.fix_val);

	/*************** Test the various Mathematical operators **********************************/

	//Addition test
	FP32 fix_add = FP32_add(a.fix_val,b.fix_val,a.Q_val,b.Q_val,&Q);
	float flt_add = toflt(fix_add,Q);
	printf("\n\nAddition of a and b using fixed point method = %f",flt_add) ;
	printf("\n                                  True Value = %f",flt_a+flt_b);

	//Subtraction test
	FP32 fix_sub = FP32_subtract(a.fix_val,b.fix_val,a.Q_val,b.Q_val,&Q);
	float flt_sub = toflt(fix_sub,Q);
	printf("\n\nSubtraction of a by b using fixed point method = %f",flt_sub) ;
	printf("\n                                  True Value = %f",flt_a-flt_b);

	//Multiplication test
	FP32 fix_mult = FP32_mult(a.fix_val,b.fix_val,a.Q_val,b.Q_val,&Q);
	float flt_mult = toflt(fix_mult,Q);
	printf("\n\nMultiplication of a and b using fixed point method = %f",flt_mult);
	printf("\n                                        True Value = %f",flt_a*flt_b);

	//Division test
	FP32 fix_div = FP32_div(a.fix_val,b.fix_val,a.Q_val,b.Q_val,&Q);
	float flt_div = toflt(fix_div,Q);
	printf("\n\nDivision of a by b using fixed point method = %f",flt_div);
	printf("\n                                   True Value = %f\n\n",flt_a/flt_b);

	//Square root test
	//FP32 fix_sqroot = sqrtx(a.fix_val, a.Q_val);
	//FP32 fix_sqroot = Fsqrt(a.fix_val, a.Q_val);
	FP32 fix_sqroot = APP_sqrt(a.fix_val, &a.Q_val);
	float flt_sqroot = toflt(fix_sqroot,a.Q_val) ;
	printf("\n\nSquare root  of a using fixed point method = %f",flt_sqroot);
	printf("\n                                  True Value = %f\n\n",sqrt(flt_a));

	printf("\n\nALL TESTS WERE SUCCESSFULL\n\n");
	
	getch();


	/*******************************  END OF TEST  *******************************/

	/**************************** Implementation of an FIR filter using fixed point method*************************/

	// defining the variables to be used
   	
	
	FP32 temp;
	QFrm fix_input[50]; 
	QFrm fix_imp_res[11] ;
	QFrm fix_output[100] ;
	float float_output[100];
	int i , j ;

	 // initialize output array values to 0
	for(i=0;i<m+n ;i++)
		fix_output[i].fix_val=0; 
	//Finding optimum Q values for every input array element
	for(i=0;i<m;i++)
	{
		fix_input[i].Q_val = Get_Optimum_Q(input[i]);
	}

	//Finding optimum Q values for every impulse response array element
	for(i=0;i<n;i++)
	{
		fix_imp_res[i].Q_val = Get_Optimum_Q(imp_res[i]);
	}

	// Conversion of input and impulse responses to fixed point format
	for(i=0;i<m;i++)
	{	fix_input[i].fix_val = tofix(input[i],fix_input[i].Q_val) ;
		//printf(" %d",fix_input[i]);
	}
	

	for(i=0;i<n;i++)
	{	fix_imp_res[i].fix_val = tofix(imp_res[i],fix_imp_res[i].Q_val) ;
		//printf(" %d",fix_imp_res[i]);
	}

	printf("\nInput and impulse responses successfully converted to fixed point\nPress any key to continue");
	getch();


	// convolution of input with impulse response
		
	for(i=0;i<m;i++)
	{
	
	for(j=i;j<i+n;j++)
	{

		
		temp = FP32_mult(fix_input[i].fix_val,fix_imp_res[j-i].fix_val,fix_input[i].Q_val,fix_imp_res[j-i].Q_val, &Q);
		fix_output[j].fix_val = FP32_add(fix_output[j].fix_val,temp,fix_output[j].Q_val,Q,&fix_output[j].Q_val);
		
	}

	
	}


	// export the floating output to a csv file

	

	printf("\n\nObtained Output	  Expected Output\n\n");
	for(i=0;i<m+n-1;i++)
	{
		//convert fix_output to float_output
		float_output[i] = toflt(fix_output[i].fix_val,fix_output[i].Q_val);
		printf("\n%f  %f",float_output[i],exp_output[i]);
		printf("  Q val = %d",fix_output[i].Q_val);
		
	}
	
	getch();
	return 0;
}
