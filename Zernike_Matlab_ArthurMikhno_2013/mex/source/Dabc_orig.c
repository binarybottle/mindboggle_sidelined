#include "mex.h"

/*
 * Dabc_orig mex file
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix 
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2006 The MathWorks, Inc.
 */

/* $Revision: 1.10.6.2 $ */

void Dabc(int facet, double *C, int N, double *D,const int *dim_C)
{
  double temp;
  int a,b,c;
  int i2,j2,k2;
  int pos1,pos2,pos3;
  int p1,p2,p3,p4,p5;
  int dims[4];
  dims[0] = *(dim_C);
  dims[1] = *(dim_C+1);
  dims[2] = *(dim_C+2);
  dims[3] = *(dim_C+3);
  
  for(a=0; a<=N; a++) {
      for(b=0; b<=N; b++) {
         for(c=0; c<=N; c++) {
             
             temp = 0;
             
             for(i2=0; i2 <= a; i2++) {
                  for(j2=0; j2 <= b; j2++) {
                       for(k2=0; k2 <= c; k2++) {
                            
                           p1 = 1;
                           p2 = i2;
                           p3 = j2;
                           p4 = k2;
                           
                           /* printf("%d , %d , % d , %d, %d \n",p1,p2,p3,p4,p5); */
                           
                           pos1 = p1+dims[0]*(p2 + dims[1]*(p3 + dims[2]*p4));
                           
                           p1 = 2;
                           p2 = a-i2;
                           p3 = b-j2;
                           p4 = c-k2;
                           
                           /* printf("%d , %d , % d , %d, %d \n",p1,p2,p3,p4,p5); */ 
                           
                           pos2 = p1+dims[0]*(p2 + dims[1]*(p3 + dims[2]*p4));
                           
                           temp = temp + *(C+pos1) * *(C+pos2);
                           
                           /* printf("%d , %d \n",pos1,pos2);            */
                           /* printf("%f , %f \n",*(C+pos1),*(C+pos2));  */
                       }
                  }
             }
            
             
             p1 = a;
             p2 = b;
             p3 = c;
             
             pos3 = p1 + (N+1)*(p2 + (N+1)*p3);
             
             *(D+pos3) = temp; 
             
         }
      }
  }
  
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *C,*D;
  double  facet;
  int N;
  int dims[3];;
  const mwSize *dim_C;
  /*int dim_array[3]; */
  
  /*  create a pointer to the input matrix y */
  C = mxGetPr(prhs[0]);
  
  /*  get the scalar order */
  N = mxGetScalar(prhs[1]);
  
  /* Get the number of dimensions in the input argument C. */
   dim_C = mxGetDimensions(prhs[0]);
  /* print ("dim_array[0] = %d", dim_array[0]); */
   
  /* Allocate the space for the return argument */
   dims[0] = N+1;
   dims[1] = N+1;
   dims[2] = N+1;
   plhs[0]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
   
  /*  create a C pointer to a copy of the output matrix */
  D = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  Dabc(facet,C,N,D,dim_C);
  
}
