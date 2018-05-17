    // Lanczos Method with Partial Reorthogonalization by S Qiao
    #include <iostream>
    #include<stdio.h>
    #include<stdlib.h>
    #include<time.h>
    #include<math.h>
	#include "Fixed_dyn.h"

//defining the definitions to be used in the program
    #define SIZ 5
    #ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif
    #define EPS 1.19209e-07
    #define SQRTEPS 3.4527e-004

// define the double type variables to be used in the program
double D1[SIZ],E1[SIZ-1],Q[SIZ][SIZ],A[SIZ][SIZ];

// 25X25
//    double A[SIZ][SIZ]={1,2,3,4,5,2,1,2,3,4,3,2,1,2,3,4,3,2,1,2,5,4,3,2,1};
	float D[SIZ],E[SIZ-1];
	typedef struct {
					FP32 fix_val;
					FP32 Q_val;
				}QFrm;

//10x10
// double A[SIZ][SIZ]={{1,2,3,4,5,6,7,8,9,2},{2,1,2,3,4,5,6,7,8,9},{3,2,1,2,3,4,5,6,7,8},{4,3,2,1,2,3,4,5,6,7},{5,4,3,2,1,2,3,4,5,6},{6,5,4,3,2,1,2,3,4,5},{7,6,5,4,3,2,1,2,3,4},{8,7,6,5,4,3,2,1,2,3},{9,8,7,6,5,4,3,2,1,2},{2,9,8,7,6,5,4,3,2,1}};
// double A[SIZ][SIZ]={1,2,3,4,5,6,7,8,9,10,2,1,2,3,4,5,6,7,8,9,3,2,1,2,3,4,5,6,7,8,4,3,2,1,2,3,4,5,6,7,5,4,3,2,1,2,3,4,5,6,6,5,4,3,2,1,2,3,4,5,7,6,5,4,3,2,1,2,3,4,8,7,6,5,4,3,2,1,2,3,9,8,7,6,5,4,3,2,1,2,10,9,8,7,6,5,4,3,2,1};

   float drand() // uniform distribution, (0..1]
    {
    return (rand()+1.0)/(RAND_MAX+1.0);
    }
    
    float random_normal() // normal distribution, centered on 0, std dev 1
    {
    return sqrt(-2*log(drand())) * cos(2*M_PI*drand())/2;
    }
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////

    void main()
    {

        srand( (unsigned)time( NULL ) );
        int i,j,k,t2,t1;
        int steps = SIZ;
        int second=0,nVec=0;
        int sqr_flag = 0 , spl_flag =0;

      //  float temp =0 ;
		QFrm temp;
       // float r[SIZ],wOld[SIZ],wCur[SIZ],qCur[SIZ],Q[SIZ][SIZ],tmp[SIZ],qOld[SIZ]; 
		QFrm r[SIZ],wOld[SIZ],wCur[SIZ],qCur[SIZ],Q[SIZ][SIZ],tmp[SIZ],qOld[SIZ];

		//Testing the optimum conversion tofix function
	/*	tofix(5.0,&temp);
		printf("\ntemp.fix_val = %d  temp.Q_val = %d",temp.fix_val,temp.Q_val);
	*/
	
        clock_t start,end; // initialize clock variable 

	/*	tofix(5.5,&r[3]);
		printf("\nr[3].fix_val = %d r[3].Q_val = %d\n", r[3].fix_val, r[3].Q_val);
		printf("r[3] actual value = %f",toflt(r[3].fix_val,r[3].Q_val));
		*/
////////////////////////////////////
 for(i=0;i<SIZ;i++)
       D1[i]=i+1;
         for(i=0;i<SIZ;i++)
		 {
			for(j=0;j<SIZ;j++)
			A[i][j]=D1[j];
			for(j=SIZ-2;j>=0;j--)
				D1[j+1]=D1[j];
				D1[0]=0;
		 }

         for(i=0;i<SIZ;i++)
		 {
         for(j=0;j<SIZ;j++)
		 {
         A[j][i]=A[i][j];
      //printf(" %f ",A[i][j]);
         }
        // printf("\n");
         }
		 printf("Matrix Formation Complete\n");
		 system("PAUSE");

		 // CONVERT THE A[][] MATRIX TO FIXED POINT FORMAT 
		 QFrm fix_A[SIZ][SIZ];
		 for(i=0;i<SIZ;i++)
		 {
         for(j=0;j<SIZ;j++)
		 {
			 fix_A[i][j].fix_val=tofix(A[i][j],fix_A[i][j].Q_val);
		 }
		 }

		 for(i=0;i<SIZ;i++){
         for(j=0;j<SIZ;j++)
			 printf("%f  ",toflt(fix_A[i][j].fix_val, fix_A[i][j].Q_val));
		 printf("\n");
		 }
		 system("PAUSE");
////////////////////////////////////
        start = clock(); 
       // Defining the r matrix
		r[0].fix_val=1;
		r[0].Q_val=14;
		wOld[0].fix_val=1;
		wOld[0].Q_val=14;
		wCur[0].fix_val=0;
		wCur[0].Q_val=14;
		qCur[0].fix_val=1;
		qCur[0].Q_val=14;
   //  initializing r matrix to N x 1 identity vector ( eye vector)
		r[0].fix_val=tofix(1.0,r[0].Q_val);
		 D[0]=0;
		 wOld[0].fix_val=tofix(1.0,wOld[0].Q_val);
		 wCur[0].fix_val=tofix(0,wCur[0].Q_val);
		 qCur[0].fix_val=tofix(1.0,qCur[0].Q_val);
		for(i=1;i<SIZ;i++)
		{ r[i].fix_val=tofix(0,r[i].Q_val);
            D[i]=0; // initialize two column vectors for diagonals
            E[i-1]=0; // Null vector formulation
			wOld[i].fix_val=tofix(0,wOld[i].Q_val);
			wCur[i].Q_val=tofix(0,wCur[i].Q_val);
			qCur[i].fix_val=tofix(0,qCur[i].Q_val);// initialize qCur
            }

		// Convert D and E matrix to their respective fixed point format
		QFrm fix_D[SIZ],fix_E[SIZ-1];

		for(i=0;i<SIZ;i++)
		{
			fix_D[i].fix_val=tofix(D[i],fix_D[i].Q_val);
			if(i<SIZ-1)
			fix_E[i].fix_val=tofix(E[i],fix_E[i].Q_val);
			//printf(" %f",toflt(fix_D[i].fix_val,fix_D[i].Q_val));
		}
// Orthogonality estimates
		FP32 t, Q_out,tt,Q_outt;
    // initialize Q
       for( i=0; i<SIZ;i++)
	   {
       Q[i][0] = qCur[i];
	//   printf("  %f",toflt(qCur[i].fix_val,qCur[i].Q_val));
	   }
    //   printf("\n");
    // Start Partial Reorthogonalization technique
    for(j=0;j<SIZ;j++)
    {// tmp = A*qCur
       for(i=0;i<SIZ;i++)
       {
		   tmp[i].fix_val=tofix(0,tmp[i].Q_val);
		   printf("Floatvalue tmp[%d]= %f",i,toflt(tmp[i].fix_val,tmp[i].Q_val));
       for(k=0;k<SIZ;k++)
		{  
			//tmp[i]+= A[i][k]*qCur[k];
			t = FP32_mult(fix_A[i][k].fix_val,qCur[k].fix_val,fix_A[i][k].Q_val,qCur[k].Q_val,&Q_out);
            printf("\nVALUE OF T %f",toflt(t,Q_out));
			tmp[i].fix_val = FP32_add(tmp[i].fix_val,t,tmp[i].Q_val,Q_out,&tmp[i].Q_val);
			printf("\ntmp[%d] = %f ",i,toflt(tmp[i].fix_val,tmp[i].Q_val));
	   }
	   
	   system("PAUSE");
       }   
	   
    // a(j) = qCur'*tmp,a[j]=0;

        for(k=0;k<SIZ;k++)
		{	
			//D[j] += qCur[k]*tmp[k];
			t = FP32_mult(tmp[k].fix_val,qCur[k].fix_val,tmp[k].Q_val,qCur[k].Q_val,&Q_out);
			fix_D[j].fix_val = FP32_add(fix_D[j].fix_val,t,fix_D[j].Q_val,Q_out,&fix_D[j].Q_val);
		//	printf("\nD[%d] = %f",j,toflt(fix_D[j].fix_val,fix_D[j].Q_val));
		}
   //   printf("\n\n\nD[%d] = %f",j,toflt(fix_D[j].fix_val,fix_D[j].Q_val));
        if(j==0)
        {
         for(k=0;k<SIZ;k++)
		 {
        // r[k] = tmp[k] - D[j]*qCur[k] ;
			 t = FP32_mult(fix_D[j].fix_val,qCur[k].fix_val,fix_D[j].Q_val,qCur[k].Q_val,&Q_out);
			 r[k].fix_val = FP32_subtract(tmp[k].fix_val,t,tmp[k].Q_val,Q_out,&r[k].Q_val);
			// printf("\nr[%d] = %f",k,toflt(r[k].fix_val,r[k].Q_val));		//correct values
		 }
         }
         else
         {
         for(k=0;k<SIZ;k++)
		 {
         //r[k] = tmp[k] - D[j]*qCur[k] - E[j-1]*qOld[k] ;
			t = FP32_mult(fix_D[j].fix_val,qCur[k].fix_val,fix_D[j].Q_val,qCur[k].Q_val,&Q_out);
			r[k].fix_val = FP32_subtract(tmp[k].fix_val,t,tmp[k].Q_val,Q_out,&r[k].Q_val);

			t = FP32_mult(fix_E[j-1].fix_val,qOld[k].fix_val,fix_E[j-1].Q_val,qOld[k].Q_val,&Q_out);
			r[k].fix_val = FP32_subtract(r[k].fix_val,t,r[k].Q_val,Q_out,&r[k].Q_val);

			//printf("\nr[%d] = %f",k,toflt(r[k].fix_val,r[k].Q_val));
		 }
         }

    if(j<SIZ-1)
    {
    //E[j] = 0;
		fix_E[j].fix_val=tofix(0,fix_E[j].Q_val);
	//	printf("\nE[%d] = %f", j,toflt(fix_E[j].fix_val,fix_E[j].Q_val) );
    for(i = 0 ; i<SIZ;i++)
	{
		//E[j] += r[i]*r[i] ;
		t = FP32_mult(r[i].fix_val,r[i].fix_val,r[i].Q_val,r[i].Q_val,&Q_out);
	//	printf("\nt = %f",toflt(t,Q_out) );
		// PROBLEM IN THIS LINE . VALUE OF E NOT GETTING UPDATED
		
		fix_E[j].fix_val = FP32_add(t,fix_E[j].fix_val,Q_out,fix_E[j].Q_val,&fix_E[j].Q_val);   
	//	printf("\nE[%d] = %f", j,toflt(fix_E[j].fix_val,fix_E[j].Q_val) );

	}

   // E[j] = sqrt(E[j]);
	fix_E[j].fix_val = APP_sqrt(fix_E[j].fix_val,&fix_E[j].Q_val);
	//printf("\n\n\nE[%d] = %f", j,toflt(fix_E[j].fix_val,fix_E[j].Q_val) );

    if(j>1)
    {

   // wOld[0]= (E[0]*wCur[1] + D[0]*wCur[0] - D[j]*wCur[0] - E[j-1]*wOld[0])/E[j] + EPS*(E[0]+E[j])*0.3*(random_normal());  
    
		 t = FP32_mult(fix_E[0].fix_val,wCur[1].fix_val,fix_E[0].Q_val,wCur[1].Q_val,&Q_out);
		 tt = FP32_mult(fix_D[0].fix_val,wCur[0].fix_val,fix_D[0].Q_val,wCur[0].Q_val,&Q_outt);

		  wOld[0].fix_val = FP32_add(t,tt,Q_out,Q_outt,&wOld[0].Q_val);

		  t = FP32_mult(fix_D[j].fix_val,wCur[0].fix_val,fix_D[j].Q_val,wCur[0].Q_val,&Q_out);

		   wOld[0].fix_val = FP32_subtract( wOld[0].fix_val,t,wOld[0].Q_val,Q_out,&wOld[0].Q_val);

		   t = FP32_mult(fix_E[j-1].fix_val,wOld[0].fix_val,fix_E[j-1].Q_val,wOld[0].Q_val,&Q_out);

			wOld[0].fix_val = FP32_subtract( wOld[0].fix_val,t,wOld[0].Q_val,Q_out,&wOld[0].Q_val);


			wOld[0].fix_val = FP32_div(wOld[0].fix_val,fix_E[j].fix_val,wOld[0].Q_val,fix_E[j].Q_val,&wOld[0].Q_val);

			// complete this part : EPS*(E[0]+E[j])*0.3*(random_normal())

	for(t2 = 1 ; t2<j-1;t2++)
    {
		//wOld[t2] = (E[t2]*wCur[t2+1] + D[t2]*wCur[t2] - D[j]*wCur[t2]+ E[t2-1]*wCur[t2-1] - E[j-1]*wOld[t2])/E[j] + EPS*0.3*(E[t2] + E[j])*random_normal();
	
		 t = FP32_mult(fix_E[t2].fix_val,wCur[t2+1].fix_val,fix_E[t2].Q_val,wCur[t2+1].Q_val,&Q_out);
		 tt = FP32_mult(fix_D[t2].fix_val,wCur[t2].fix_val,fix_D[t2].Q_val,wCur[t2].Q_val,&Q_outt);

		 wOld[t2].fix_val = FP32_add(t,tt,Q_out,Q_outt,&wOld[t2].Q_val);

		 t = FP32_mult(fix_D[j].fix_val,wCur[t2].fix_val,fix_D[j].Q_val,wCur[t2].Q_val,&Q_out);

		 wOld[t2].fix_val = FP32_subtract( wOld[t2].fix_val,t,wOld[t2].Q_val,Q_out,&wOld[t2].Q_val);
		 
		 t = FP32_mult(fix_E[t2-1].fix_val,wCur[t2-1].fix_val,fix_E[t2-1].Q_val,wCur[t2-1].Q_val,&Q_out);

		 wOld[t2].fix_val = FP32_add( wOld[t2].fix_val,t,wOld[t2].Q_val,Q_out,&wOld[t2].Q_val);

		 t = FP32_mult(fix_E[j-1].fix_val,wOld[t2].fix_val,fix_E[j-1].Q_val,wOld[t2].Q_val,&Q_out);

		 wOld[t2].fix_val = FP32_subtract( wOld[t2].fix_val,t,wOld[t2].Q_val,Q_out,&wOld[t2].Q_val);

		 wOld[t2].fix_val = FP32_div(wOld[t2].fix_val,fix_E[j].fix_val,wOld[t2].Q_val,fix_E[j].Q_val,&wOld[t2].Q_val);

		 // complete this part : EPS*(E[t2]+E[j])*0.3*(random_normal())
	}
	
	//printf("\nwOld[5] = %d",wOld[5].fix_val);
	for(t2 = 0 ; t2< j-1 ;t2++)
    {
    temp = wOld[t2] ;
    wOld[t2] = wCur[t2] ;
    wCur[t2] = temp ;
	wOld[j-1].fix_val=tofix(1.0,wOld[j-1].Q_val); 
    }
    }// if(j>1)
	//wCur[j] = EPS*SIZ*(E[1]/E[j])*0.6*(random_normal())
	wCur[j].fix_val=tofix(EPS*SIZ*((toflt(fix_E[1].fix_val,fix_E[1].Q_val)/toflt(fix_E[j].fix_val,fix_E[j].Q_val)))*0.6*(random_normal()),wCur[j].Q_val ) ;
    //wCur[j+1] = 1.0 ;
	wCur[j+1].fix_val=tofix(1.0,wCur[j+1].Q_val);
    
    for(t2=0;t2<j;t2++)
    {
		if( (toflt(wCur[t2].fix_val,wCur[t2].Q_val) > SQRTEPS) || (toflt(wCur[t2].fix_val,wCur[t2].Q_val) < - SQRTEPS ))
    {
    sqr_flag = 1;
    spl_flag =1;
    }
    }
    if( second == 1 )
    {
    sqr_flag = 1;
    }
    if(spl_flag ==1 )
    second = 1;
    else
    second = 0;
    if(sqr_flag ==1)
    {
    for(k=0;k<j;k++)
    {
		temp.fix_val=tofix(0,temp.Q_val);
	 for(t1=0;t1<SIZ;t1++)
	 {
		// temp += Q[t1][k]*r[t1];
			 t = FP32_mult(Q[t1][k].fix_val,r[t1].fix_val,Q[t1][k].Q_val,r[t1].Q_val,&Q_out);
			 temp.fix_val = FP32_add(temp.fix_val,t,temp.Q_val,Q_out,&temp.Q_val);
	 }
    for(t2=0;t2<SIZ;t2++)
    {
     //r[t2] = r[t2] - temp*Q[t2][k];
		
		t = FP32_mult(Q[t2][k].fix_val,temp.fix_val,Q[t2][k].Q_val,temp.Q_val,&Q_out);
		r[t2].fix_val = FP32_subtract(r[t2].fix_val,t,r[t2].Q_val,Q_out,&r[t2].Q_val);

    }
     //wCur[k] = EPS*1.5*(random_normal());
	wCur[k].fix_val=tofix(EPS*1.5*(random_normal()),wCur[k].Q_val);
    }
    nVec = nVec + j ;
	fix_E[j].fix_val=tofix(0,fix_E[j].Q_val);

    for(i = 0 ; i<SIZ; i++)
	{ 
		//E[j] += r[i]*r[i] ;
		t = FP32_mult(r[i].fix_val,r[i].fix_val,r[i].Q_val,r[i].Q_val,&Q_out);
		fix_E[j].fix_val = FP32_add(fix_E[j].fix_val,t,fix_E[j].Q_val,Q_out,&fix_E[j].Q_val);
	}

   // E[j] = sqrt(E[j]);
	fix_E[j].fix_val = APP_sqrt(fix_E[j].fix_val,&fix_E[j].Q_val);
	

    }// if(sqr_flag ==1)
    if(j+1 < EPS)
    {
    steps = j ;
    break;
    }
    else
    {
    for(t2 = 0;t2<SIZ ; t2++)
    {  
		//qOld[t2] = qCur[t2];
    //qCur[t2] = r[t2]/E[j] ;
   // Q[t2][j+1] = qCur[t2];
    qOld[t2] = qCur[t2];
    qCur[t2].fix_val = FP32_div(r[t2].fix_val,fix_E[j].fix_val,r[t2].Q_val,fix_E[j].Q_val,&qCur[t2].Q_val) ;
    Q[t2][j+1] = qCur[t2];
    }
	
    }
    }
    }
    // TT = form the tridiagonalized matrix ///////////////////////////////// D n C ///
	
//////////////////////////////////////////////////////////
//	printf("\nD[5] = %f",toflt(fix_D[5].fix_val,fix_D[5].Q_val));
    end = clock();
	printf("\n\nTotal elapsed time : %d ms\n",(end-start));
    printf("\nPrimary Diagonals Elements::\n\n");
    for(i=0;i<SIZ;i++)
		printf("%f,",toflt(fix_D[i].fix_val,fix_D[i].Q_val));
    printf("\n");
    printf("\nOff Diagonal Elements::\n\n");
    for(i=0; i<SIZ-1;i++)
		printf("%f,",toflt(fix_E[i].fix_val,fix_E[i].Q_val));
  //  printf("\n Q Orthogonal Matrix");
//    printf("\n");
 /*  for(i=0;i<SIZ;i++)
    {
                    for(j=0;j<SIZ;j++)
                    printf(" %f ",toflt(Q[i][j].fix_val,Q[i][j].Q_val) );
                    printf("\n");
} */
	
    //printf("\nstart = %d end = %d",start,end);
   
	 
	system("PAUSE");
    
	}
	
	
///////////////////////////////////////////////////////////////////////////////
	