typedef unsigned long int FP32 ;


FP32 fconv( FP32 a , FP32 q1 , FP32 q2);
int faddg(int a , int b , int q1 , int q2 , int q3);
FP32 tofix(float a, FP32 q);
float toflt(FP32 a , FP32 q);


#define	INF_FP32	0x7fffffffUL
#define	NINF_FP32	0xffffffffUL

#define	FP32_ZERO	0x00000000UL
#define FP32_MSB	0x80000000UL

#define	FP32_LO_WORD(x)		(*((unsigned short int *) (&x)))		// the type of FP32_32 is ULL, so it divides integer part to HI word,
#define	FP32_HI_WORD(x)		(*((unsigned short int *) (&x) + 1))	// and exponetial part to LO word. that why he used UL to represent value.
#define	FP32_SIGN(x)			(x & 0x80000000UL)	
#define	FP32_CLEAR_SIGN(x)	(x &= 0x7fffffffUL)	// 
#define	FP32_SET_SIGN(x)		(x|=0x80000000UL)  // 
#define	FP32_PUT_SIGN(x,y)	if(y){FP32_SET_SIGN(x);}else{FP32_CLEAR_SIGN(x);}
#define	FP32_TOGGLE_SIGN(x)	(x^=0x80000000)
#define	FP32_INTEGER(x)		(x & 0x7fffffffUL)>>Q
#define	ONE		0x00000001UL

FP32 FP32_add (FP32 num1, FP32 num2, FP32 Q1, FP32 Q2, FP32 *Q) ;

FP32 FP32_subtract (FP32 num1, FP32 num2, FP32 Q1, FP32 Q2, FP32 *Q) ;
FP32 FP32_mult (FP32 num1, FP32 num2, FP32 Q1, FP32 Q2, FP32 *Q) ;
FP32 FP32_div (FP32 num1, FP32 num2, FP32 Q1, FP32 Q2, FP32 *Q) ;
FP32 sqrtx (FP32 a, FP32 Q) ;
FP32 Get_Optimum_Q(float flt_a);
FP32 Fsqrt(FP32 M , FP32 Q);
FP32 adjust(FP32 x, FP32 Q);
FP32 APP_sqrt(FP32 x , FP32* Q);