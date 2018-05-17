#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "Fixed_dyn.h"

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

float imp_res[11] = {0.061036,0.074909,0.086908,0.096172,0.10202,0.10402,0.10202,0.096172,0.086908,0.074909,0.061036};

//Output for Input (1) obtained from Matlab
float exp_output[61] = {0,0.014825,0.044983,0.088065,0.13991,0.19628,0.25446,0.31401,0.37656,0.4446,0.51984,0.60171,0.67545,0.74073,0.79653,0.84112,0.87245,0.88866,0.88866,0.87245,0.84112,0.79653,0.74073,0.67545,0.60171,0.51984,0.42978,0.33157,0.22594,0.11455,2.498e-016,-0.11455,-0.22594,-0.33157,-0.42978,-0.51984,-0.60171,-0.67545,-0.74073,-0.79653,-0.84112,-0.87245,-0.88866,-0.88866,-0.87245,-0.84112,-0.79653,-0.74073,-0.67545,-0.60171,-0.51984,-0.4446,-0.37656,-0.31401,-0.25446,-0.19628,-0.13991,-0.088065,-0.044983,-0.014825,-2.9899e-017};
int main()
{
	
	FP32 Q;// = 14 ;
	/**************************** Implementation of an FIR filter using fixed point method*************************/

//	 defining the variables to be used
   	
	
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
	
	system("pause");
	return 0;
}