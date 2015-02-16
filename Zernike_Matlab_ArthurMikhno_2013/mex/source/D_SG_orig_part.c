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

void D_SG_orig_part(int num_facets, int i, int j, int k, double *C, double *D, double *S, double *Vol, double *F,
        const int *dim_C, const int *dim_D, const int *dim_S)
{
    int i1,j1,k1,facet;
    int pos1,pos2,pos3;
    int p1,p2,p3,p4,p5;
    int dimC[4], dimD[4], dimS[4];
    double aux_1, aux_2, aux, tmp;
    
    dimC[0] = *(dim_C);
    dimC[1] = *(dim_C+1);
    dimC[2] = *(dim_C+2);
    dimC[3] = *(dim_C+3);
    
    dimD[0] = *(dim_D);
    dimD[1] = *(dim_D+1);
    dimD[2] = *(dim_D+2);
    dimD[3] = *(dim_D+3);
    
    dimS[0] = *(dim_S);
    dimS[1] = *(dim_S+1);
    dimS[2] = *(dim_S+2);
    dimS[3] = *(dim_S+3);
    
    
    /* pre-compute factorials */
    
    
    for (facet=0; facet < num_facets; facet++) {
        aux_1 = *(F + i) * *(F + j) * *(F + k);
        aux_2 = *(F + i+j+k+2);
        aux = aux_1 / aux_2;
        
        tmp = 0;
        
        for(i1=0; i1 <= i; i1++) {
            for(j1=0; j1 <= j; j1++) {
                for(k1=0; k1 <= k; k1++) {
                    
                    p1 = facet;
                    p2 = i1;
                    p3 = j1;
                    p4 = k1;
                    
                    /* printf("%d , %d , % d , %d, %d \n",p1,p2,p3,p4,p5); */
                    
                    pos1 = p1+dimC[0]*(p2 + dimC[1]*(p3 + dimC[2]*p4));
                    
                    p1 = facet;
                    p2 = i-i1;
                    p3 = j-j1;
                    p4 = k-k1;
                    /* printf("%d , %d , % d , %d, %d \n",p1,p2,p3,p4,p5); */
                    
                    pos2 = p1+dimD[0]*(p2 + dimD[1]*(p3 + dimD[2]*p4));
                    
                    tmp = tmp + *(C + pos1) * *(D + pos2);
                    
                }
            }
        }
        
        tmp = tmp * aux;
        
        p1 = facet;
        p2 = i;
        p3 = j;
        p4 = k;
        /* printf("%d , %d , % d , %d, %d \n",p1,p2,p3,p4,p5); */
        pos3 = p1+dimS[0]*(p2 + dimS[1]*(p3 + dimS[2]*p4));
        *(S + pos3) = (tmp * *(Vol + facet)) / (i+j+k+3);
        
    }
    
    
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *C,*D,*Vol,*S,*F;
    int num_facets,i,j,k;
    int dims[3];;
    const mwSize *dim_C;
    const mwSize *dim_D;
    const mwSize *dim_S;
    
    
    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgTxt
     * within an if statement, because it will never get to the else
     * statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     * the MEX-file) */
    
    
    /*  get the scalar input x */
    num_facets = mxGetScalar(prhs[0]);
    i          = mxGetScalar(prhs[1]);
    j          = mxGetScalar(prhs[2]);
    k          = mxGetScalar(prhs[3]);
    
    /*  create a pointer to the input matrix y */
    C = mxGetPr(prhs[4]);
    D = mxGetPr(prhs[5]);
    /*Stemp = mxGetPr(prhs[6]); */
    Vol = mxGetPr(prhs[7]);
    F   = mxGetPr(prhs[8]);
    
    
    /* Get the number of dimensions in the input argument C. */
    dim_C = mxGetDimensions(prhs[4]);
    dim_D = mxGetDimensions(prhs[5]);
    dim_S = mxGetDimensions(prhs[6]);
    
    
    /* Allocate the space for the return argument */
    dims[0] = *(dim_S);
    dims[1] = *(dim_S+1);
    dims[2] = *(dim_S+2);
    dims[3] = *(dim_S+3);
    plhs[0]=mxCreateNumericArray(4,dims,mxDOUBLE_CLASS,mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    S = mxGetPr(plhs[0]);
    
    /*  call the C subroutine */
    D_SG_orig_part(num_facets,i,j,k,C,D,S,Vol,F,dim_C,dim_D,dim_S);
    
}
