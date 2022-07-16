/*=================================================================
 * The main routine analyzes all incoming (right-hand side) arguments
 *
 * Copyright 1984-2011 The MathWorks, Inc.
 *
 *=================================================================*/

// #include <stdio.h>
// #include <string.h>
#include "mex.h"
// input: UniqueG,IndexQ,IndexP,Indexf,n
long long CZAdd(long long x,long long y ,long long n2){
    return n2-((n2-x)^(n2-y));
}
void mexFunction(int nlhs, mxArray *plhs[] ,int nrhs, const mxArray *prhs[])//(int nlhs, mxArray *plhs[])
{
    double *Indexf,*IndexQ,*IndexP,*fhat,*nn,*UniqueG;
    long long n;
    long long *indf, *indg;
    nn=mxGetPr(prhs[3]);
    n=(long long)nn[0];
    IndexQ=mxGetPr(prhs[0]);
    IndexP=mxGetPr(prhs[1]);
    Indexf=mxGetPr(prhs[2]);
    int lengthQ=(int)mxGetM(prhs[0]);
    int lengthP=(int)mxGetM(prhs[1]);// number of rows  of var lengthP
    int lengthf=(int)mxGetM(prhs[2]);
    long long n2=1;
    for (int i=0;i<n;i++){
        n2=2*n2;
    }
    long long tempn;
    long long *FProdQ= ( long long*)mxMalloc( sizeof(long long)*lengthf*lengthQ*lengthQ*4) ;
    long long *indP= ( long long*)mxMalloc( sizeof(long long) * lengthP*lengthP );
    for(int i=0;i<lengthf;i++){
        for(int j=0;j<lengthQ;j++){
            for(int k=0;k<lengthQ;k++){
                // Eq1[w*lengthQ*lengthQ+j*lengthQ+k]+=fhat[i];// (w*lengthQ*lengthQ+i*lengthQ+j->(w,i,j))
                FProdQ[i*lengthQ*lengthQ+j*lengthQ+k]=CZAdd(CZAdd((long long)IndexQ[k],(long long)IndexQ[j],n2),(long long)Indexf[i],n2);
            }
        }
    }
    for(int i=0;i<lengthP;i++){
        for(int j=0;j<lengthP;j++){
            indP[i+j*lengthP]=CZAdd((long long)IndexP[i],(long long)IndexP[j],n2);
        }
    }
    plhs[0] = mxCreateDoubleMatrix(lengthf*lengthQ*lengthQ ,4, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(lengthP*lengthP ,3, mxREAL);
    double *output1,*output2;
    output1=(double*)mxGetData(plhs[0]);
    output2=(double*)mxGetData(plhs[1]);
    mexPrintf("\n");
    int LengthG=lengthf*lengthQ*lengthQ;
    for(int i=0;i<lengthf;i++){
        for(int j=0;j<lengthQ;j++){
            for(int k=0;k<lengthQ;k++){
                output1[i*lengthQ*lengthQ+j*lengthQ+k]=(double)FProdQ[i*lengthQ*lengthQ+j*lengthQ+k];
                output1[(i*lengthQ*lengthQ+j*lengthQ+k)+LengthG]=(double)i+1;
                output1[(i*lengthQ*lengthQ+j*lengthQ+k)+LengthG*2]=(double)j+1;
                output1[(i*lengthQ*lengthQ+j*lengthQ+k)+LengthG*3]=(double)k+1;
            }
        }
    }
    LengthG=lengthP*lengthP;
    for(int i=0;i<lengthP;i++){
        for(int j=0;j<lengthP;j++){
            output2[i+j*lengthP+LengthG]=(double)i+1;
            output2[i+j*lengthP+LengthG*2]=(double)j+1;
           output2[i+j*lengthP]=indP[i+j*lengthP];
        }
    }
    mxFree(FProdQ);
    mxFree(indP);
    return;
}