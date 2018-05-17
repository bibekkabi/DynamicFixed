/***************************************************************************
Various mathematical operations defined for fixed point ALU

***************************************************************************/

#include<stdio.h>

/* The basic operations perfomed on two numbers a and b of fixed
point q format returning the answer in q format */
#include <tchar.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <conio.h>

#include "Fixed_point.h"


unsigned int g_uiOverflowAddition = 0 ;								// Addition overflow
unsigned int g_uiOverflowMultiplication = 0 ;						// Multiplication overflow
unsigned int g_uiOverflowDivision = 0 ;								// Division overflow
//FP32 Q;

FP32 FP32_add (FP32 num1, FP32 num2, FP32 Q1, FP32 Q2,FP32 *Q)
{    
	FP32 Qtmp;
	//printf("\nQ = %d",*Q);
	//getch();
	if(Q1<Q2)
	{
		*Q = Q1-1;
		
	}
	else 
	{
		*Q = Q2-1;
	}
    /*printf("\nBefor fix_a = %f   fix_b = %f\n",toflt(num1,Q1),toflt(num2,Q2));
	printf("\nQ1=%d,Q2=%d and Q=%d",Q1,Q2,*Q);*/
	num1 = fconv(num1,Q1,*Q);
	num2 = fconv(num2,Q2,*Q);
	/*printf("\nModified Q = %d",*Q);
	printf("\n After fix_a = %ld   fix_b = %ld\n",num1,num2);*/
	/*getch();*/
	FP32 result_sign ;
	FP32 sum ;

	if (FP32_SIGN(num1) == FP32_SIGN(num2))
	{
		// perform addition (num1 + num2)
		result_sign = FP32_SIGN(num1) ;
		sum = num1 + num2 ;
		
		/*
			num1 and num2 both are positive or both are negative
			i.e. their MSBits are both 0 or both 1.
			We want to check if a carry is generated from bit62
			i.e. the bit next to MSBit.
			If carry is NOT generated:
			then MSBit of sum will be 0 (0 + 0 = 0 OR 1 + 1 = 0)
			If carry is generated:
			then MSBit of sum will be 1 (0 + 0 + 1 = 1 OR 1 + 1 + 1 = 1)
			So if MSBit of sum is set, it means there is overflow!
		*/
		if (FP32_SIGN(sum))
		{
			// overflow error
			g_uiOverflowAddition = 1 ;
			printf("\n\n !!!!OVERFLOW!!!!\n\n");
			getch();
		}
	}
	else
	{
		// perform subtraction (num1 - num2)
		// before subtracting, we convert the negative number to positive
		// i.e. ultimately, always subtraction is performed on 2 positive numbers
		// besides, always (large_number - small_number) is performed
		// so there will never be any borrow
		// i.e. no need to check for overflow
		if (FP32_SIGN(num2))
		{
			// num1 positive, num2 negative
			//printf("b negative"); 
			FP32_CLEAR_SIGN(num2) ;
			//printf("\n\nb = %d\n",num2);
			if (num1 >= num2)
			{
				// absolute value of num1 is >= that of num2
				// hence, result will be positive
				result_sign = 0 ;
				sum = num1 - num2 ;
			
			}
			else
			{
				// absolute value of num1 is < that of num2
				// hence, result will be negative
				result_sign = 1 ;
				sum = num2 - num1 ;
			}
		}
		else
		{
			// num1 negative, num2 positive
			FP32_CLEAR_SIGN(num1) ;
			if (num2 >= num1)
			{
				// absolute value of num2 is >= that of num1
				// hence, result will be positive
				result_sign = 0 ;
				sum = num2 - num1 ;
			}
			else
			{
				// absolute value of num2 is < that of num1
				// hence, result will be negative
				result_sign = 1 ;
				sum = num1 - num2 ;
			}
		}
	}
	FP32_PUT_SIGN(sum,result_sign) ;
	Qtmp=Get_Optimum_Q(toflt(sum,*Q));
	sum=tofix(toflt(sum,*Q),Qtmp);
	*Q=Qtmp;
	//printf("\nresult=%f\n and Optimum Q=%d and returned Q=%d",toflt(sum,*Q),Qtmp,*Q);
	//getch();
	return sum ;
}

FP32 FP32_subtract (FP32 num1, FP32 num2 , FP32 Q1, FP32 Q2, FP32 *Q)
{
	FP32_TOGGLE_SIGN(num2) ;
	return FP32_add(num1,num2,Q1,Q2,Q) ;
}


FP32 FP32_mult (FP32 num1, FP32 num2, FP32 Q1, FP32 Q2, FP32 *Q)
{
	FP32 product ;
	FP32 ullTemp ;
	FP32 result_sign ;
	FP32 Qtmp2;

	if((31-Q1+31-Q2) < 31)
	{
		*Q = 31-(31-Q1+31-Q2);
	}
	else
		*Q=0;
    /*printf("\n Before fix_a = %f   fix_b = %f\n",toflt(num1,Q1),toflt(num2,Q2));
	printf("\n Q1=%ld and Q2=%ld\n",Q1,Q2);
	printf("\nBefore Modification Q = %d",*Q);*/
	num1 = fconv(num1,Q1,*Q);
	num2 = fconv(num2,Q2,*Q);
    /*printf("\nModified Q = %d",*Q);
	printf("\n After fix_a = %ld   fix_b = %ld\n",num1,num2);
	getch();*/
	
	result_sign = FP32_SIGN(num1) ^ FP32_SIGN(num2) ; // XOR operation to find the resulatant sign
	FP32_CLEAR_SIGN(num1) ;
	FP32_CLEAR_SIGN(num2) ;


	// in multiplication , we separate the number into lo_word and hi_word, then it doesn't need to cosider more bit place for our number
	// moreover, while it calculates different Q format, it only shift the difference of binary point with Q16
	product = ((unsigned long int) FP32_LO_WORD(num2)) * ((unsigned long int) FP32_HI_WORD(num1)) ;
	if(*Q>=16)
	product >>= (*Q-16);
	else
	product <<= (16-*Q);

	if (product & 0x80000000UL)		//overflow
	{
		g_uiOverflowMultiplication = 1 ;
		//printf("!!!OVERFLOW!!!");
		//getch();
	}


	ullTemp = ((unsigned long int) FP32_HI_WORD(num2)) * ((unsigned long int) FP32_LO_WORD(num1)) ;
	if(*Q>=16)
	ullTemp >>= (*Q-16);
	else
	ullTemp <<= (16-*Q);
	
	if (ullTemp & 0x80000000UL)
	{
		g_uiOverflowMultiplication = 1 ;
		//printf("!!!OVERFLOW!!!");
		//getch();
	}
	
	product += ullTemp ;
	if (product & 0x80000000UL)
	{
		g_uiOverflowMultiplication = 1 ;
		//printf("!!!OVERFLOW!!!");
		//getch();
	
	}

	// While two fractional number are mutiplied, it should be sfited back Q bit to get correct frational position.
	ullTemp = ((unsigned long int) FP32_LO_WORD(num2)) * ((unsigned long int) FP32_LO_WORD(num1)) ;
	product += ullTemp >> *Q;
	if (product & 0x80000000UL)
	{
		g_uiOverflowMultiplication = 1 ;
		//printf("!!!OVERFLOW!!!");
		//getch();
	
	}
	if (ullTemp & 0x80000000UL)
	{
		product ++ ;
		if (product & 0x80000000UL)
		{
			g_uiOverflowMultiplication = 1 ;
			//printf("!!!OVERFLOW!!!");
			//getch();
	
		}
	}

	//While two integer numbers are multiplied , the product should be left shift to the correct integer position/.
	ullTemp = ((unsigned long int) FP32_HI_WORD(num2)) * ((unsigned long int) FP32_HI_WORD(num1)) ;
	if (ullTemp & 0xffff8000UL)
	{
		g_uiOverflowMultiplication = 1 ;
		//printf("!!!OVERFLOW!!!");
		//getch();
	
	}
	product += (ullTemp << (32-*Q) ) ;

	if (product & 0x80000000UL)
	{
		g_uiOverflowMultiplication = 1 ;
		//printf("!!!OVERFLOW!!!");
		//getch();
	}

	// Finally, u have to put the result sign back to the profuct.
	FP32_PUT_SIGN(product,result_sign) ;
    Qtmp2=Get_Optimum_Q(toflt(product,*Q));
	product=tofix(toflt(product,*Q),Qtmp2);
	*Q=Qtmp2;
	return product ;

}



FP32 FP32_div (FP32 num1, FP32 num2, FP32 Q1, FP32 Q2, FP32 *Q)
{
	FP32 result_sign, ullTemp , q ;
	FP32 Qtmp3;
	unsigned short int u ;

	//Condition to choose the optimum Q value for division of num1 by num2
	if(((num1>>Q1) >=1) && ((num2>>Q2)>=1))
	{
		if(Q1>Q2)
			*Q = Q2;
		else
			*Q = Q1;
	}
	else if((num2>>Q2) <1)
	{
		if((Q1+Q2)<31)
			*Q = Q1+Q2;
		else 
			*Q = 16;
	}
	num1 = fconv(num1,Q1,*Q);
	num2 = fconv(num2,Q2,*Q);

	result_sign = FP32_SIGN(num1) ^ FP32_SIGN(num2) ;
/*
	numerator == 0
	denominator == 0
	q too large
	q too small
*/
	FP32_CLEAR_SIGN(num1) ;
	FP32_CLEAR_SIGN(num2) ;
	if ((num2 == 0) || ((ullTemp = num1/num2) > 0x7fffffffUL))  // At here, u already get the integer value of the quotient
	{
		g_uiOverflowDivision = 1 ;
		if (result_sign)
		{
			return NINF_FP32 ;
		}
		else
		{
			return INF_FP32 ;
		}
	}
	q = ullTemp << *Q ;	// we have found the integer part of quotient already

	ullTemp = num1 % num2 ;	// ulTemp will be < num2, but still maybe more than 32 bits, so here it has to calculate the quotient of fractional part
	for(u = 0 ; u < *Q ; u ++)	// here, u have to shift the fractional part of number untill 0x80000000UL , then divid it.
	{
		if (ullTemp & 0x80000000UL)
		{
			break ;
		}
		ullTemp <<= 1 ;	
	}
	if (u < *Q)
	{
		num2 >>= (*Q - u) ;
	}
	
	q += (ullTemp / num2) ;				// fractional part of quotient
	FP32_PUT_SIGN(q,result_sign) ;
    Qtmp3=Get_Optimum_Q(toflt(q,*Q))-2;
	q=tofix(toflt(q,*Q),Qtmp3);
	*Q=Qtmp3;
	return q ;
}

FP32 sqrtx (FP32 a  , FP32 Q)
{
	
// this algorithm is download from the website http://www.codecodex.com/wiki/Calculate_an_integer_square_root
	// the adbvantage of this algorithm that we don't need to use multiplication and division, 
	// it only uses addition and scaling, so the overflow problem is not as serious as above Newton's method.

	register FP32 // OR register uint16 OR register uint8 - respectively  
     root, remainder, place;  
  
    root = 0;  
    remainder = a;  
    place = 0x40000000; // OR place = 0x4000; OR place = 0x40; - respectively  
  
    while (place > remainder)  
        place = place >> 2;  
    while (place)  
    {  
        if (remainder >= root + place)  
        {  
            remainder = remainder - root - place;  
            root = root + (place << 1);  
        }  
        root = root >> 1;  
        place = place >> 2;  
    }  
    return root<<(Q/2);  //integer square root shift ur Q/2 back
}



FP32 Fsqrt(FP32 M , FP32 Q) //square-root calculation of Fixed datatype
{
if(Q%2==1)
{(M = M>>1);Q--;}
register FP32 // OR register uint16 OR register uint8 - respectively
root, remainder, place;

root = 0;
remainder = M;
place = 0x40000000; // OR place = 0x4000; OR place = 0x40; - respectively

while (place > remainder)
place = place >> 2;
while (place)
{
if (remainder >= root + place)
{
remainder = remainder - root - place;
root = root + (place << 1);
}
root = root >> 1;
place = place >> 2;
}
M=root<<(Q/2);
return adjust(M,Q/2);
}

FP32 adjust(FP32 x , FP32 Q)
{ FP32 integer, n=ONE; 
int i=0;
integer=x;
if (FP32_SIGN(integer)){
FP32_CLEAR_SIGN(integer);
}
if(integer==0) Q=31;
else
{
while(integer>=n)
{++i;
n=n<<1; }
if(i>=Q)
x=fconv(x,Q,30-i+Q);
else
x=fconv(x,Q,30);
}
return x;
}


/* convert a from q1 format to q2 format */

FP32 fconv( FP32 a , FP32 q1 , FP32 q2)
{
	if(FP32_SIGN(a))
	{
		//printf("\nsign change\n");
		FP32_TOGGLE_SIGN(a);

		if(q2>q1)
		{
			a = a<<q2-q1;
		}
		else
		{
			a = a>>q1-q2;
		}
		//((q2)>(q1)) ? (a)<<((q2)-(q1)) : (a)>>((q1)-(q2)) ;
		FP32_TOGGLE_SIGN(a);
		return a;

	}
	else
	{
	//	printf("\nNo sign change\n");
		if(q2>q1)
		{
			a = a<<q2-q1;
		}
		else
		{
			a = a>>q1-q2;
		}
		//return(((q2)>(q1)) ? (a)<<((q2)-(q1)) : (a)>>((q1)-(q2))) ;
	}
}

/* the general operation between a in q1 format and b in q2 format
returning the result in q3 format */

int faddg(int a , int b , int q1 , int q2 , int q3)
{
	return ( fconv(a,q1,q3) + fconv(b,q2,q3) ) ;
}

/* convert to and from floating point */

FP32 tofix(float a, FP32 q)
{
	
	if(a>=0)
	return ((FP32)( (a)*(float)(1<<(q)) )) ;
	else 
	{
		a= -a;
		FP32 x = ((unsigned int)( (a)*(float)(1<<(q)) )) ;
		FP32_SET_SIGN(x);
		return x ;
		
	}

}

float toflt(FP32 a ,FP32 q) 
{
	int x ;
	if(FP32_SIGN(a))
	{
		FP32_CLEAR_SIGN(a);
		x = a;
		x = -x;
		
	}
	else
		x=a;
	return ( (float)(x) / (float)(1<<(q)) ) ;
}

FP32 Get_Optimum_Q(float flt_a)
{
	FP32 q_val =1;
	FP32 Qtest = 2;
	int tshift;
	tshift=((unsigned int(abs(flt_a)) != 0) && ((unsigned int(abs(flt_a)) & (~unsigned int(abs(flt_a)) + 1)) == unsigned int(abs(flt_a))));
	
	if(tshift==1){
		
        tshift= log10(abs(flt_a))/log10(2.0);
		q_val=q_val+tshift;
		//printf("\nMultiples of 2\n");
	}
	else{
	while(abs(flt_a)>Qtest)
	{		Qtest = Qtest<<1;
			q_val++;	
	}
	//printf("\nOther numbers\n");
	}
	//printf("value=%f and q_val = %d\n",flt_a,q_val);
	return (31 - q_val);
}




//APPLE SQUARE ROOT ALGORITHM

FP32 APP_sqrt(FP32 x , FP32 *Q)
{
 unsigned long root, remHi, remLo, testDiv, count;
root = 0; /* Clear root */
remHi = 0; /* Clear high part of partial remainder */
remLo = x; /* Get argument into low part of partial remainder */
count = 30; /* Load loop counter */
FP32 shift = 0;
int Qtmp1;
if(FP32_SIGN(x))
{
	printf("\nSQUARE ROOT OF A NEGATIVE NUMBER BEING CALCULATED !!");
	getch();
	return NULL;
}

//FP32 temp =x ;
/*if(*Q%2 !=0)
{
	x = x>>1;
	*Q = *Q-1;
}
*/
x = x>>(30-*Q);
/*if(*Q%2 ==0)
shift = (30-*Q)/2;
else
shift = (29-*Q)/2;
*/
//shift = (*Q%2==0)?(30-*Q)/2:(30-*Q-1)/2;
remLo = x;
//printf("\nshift = %d",shift);
do {
remHi = (remHi<<2) | (remLo>>30); remLo <<= 2; /* get 2 bits of arg */
root <<= 1; /* Get ready for the next bit in the root */
testDiv = (root << 1) + 1; /* Test radical */
if (remHi >= testDiv) {
remHi -= testDiv;
root++;
}
} while (count-- != 0);
//printf("Root result =%f and Q=%d",toflt)
//Qtmp1=Get_Optimum_Q(toflt(root,*Q))-1;
//root=tofix(toflt(root,*Q),Qtmp1);
//*Q=Qtmp1;
return(root);
}