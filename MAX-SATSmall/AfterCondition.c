//function After
// // // // // Input:
// SE1,SE2 was sorted
// UniqueG is also sorted
// [subs,val]=find(f)
// Np=length(P); Nq=length(Q);
// // // // // Output: Index Vector IndexSP,IndexSQ,PKStart,PKEnd,QKStart,QKEnd and double array F s.t. the condition of SDP is
// // min(norm(T,1))
// // for i=1:length(UniqueG)
// //       T(i)=dot(F(QKStart(i):QKEnd(i)),Q(IndexSQ(QKStart(i):QKEnd(i))))-sum([0,P(IndexSP(PKStart(i):PKEnd(i)))]);
// // end
#include "mex.h"
double Fmap(double x, double *subs , double *val, long long nsubs){
    double z=0;
    for(int i=0;i<nsubs;i++){
        if(subs[i]==x){
            z=val[i];
            break;
        }
    }
    return z;
    
}
int min(int x,int y ){
    if (x<y){
        return x;
    }
    else{
        return y;
    }
}
void mexFunction(int nlhs, mxArray *plhs[] ,int nrhs, const mxArray *prhs[])//(int nlhs, mxArray *plhs[])
{
    //[IndexSP,IndexSQ,F,PKStart,PKEnd,QKStart,QKEnd]=after(SE1,SE2,UniqueG,subs,val,Np,Nq)
    double *SE1, *SE2,*doubleUniqueG, *doublesubs,*val,*doubleNp,*doubleNq;
    // double *Indexf,*IndexQ,*IndexP,*fhat,*nn,*UniqueG;
    int LengthG=mxGetM(prhs[2])*mxGetN(prhs[2]);
    int Lengthf=mxGetM(prhs[4])*mxGetN(prhs[4]);
    //  long long *UniqueG= ( long long*)mxMalloc( sizeof(long long)*LengthG) ;
    // long long *subs=( long long*)mxMalloc( sizeof(long long)*Lengthf) ;
    int Np, Nq;
    //取出输入参数
    SE1=mxGetPr(prhs[0]);
    SE2=mxGetPr(prhs[1]);
    
    doubleUniqueG=mxGetPr(prhs[2]);
    doublesubs=mxGetPr(prhs[3]);
    val=mxGetPr(prhs[4]);
    doubleNp=mxGetPr(prhs[5]);
    doubleNq=mxGetPr(prhs[6]);
    Np=(int)doubleNp[0];
    Nq=(int)doubleNq[0];
    double *F=(double *)mxMalloc( sizeof(double)*mxGetM(prhs[0])) ;
    int *QKStart=(int*)mxMalloc( sizeof(int)*LengthG) ;
    int *QKEnd=(int*)mxMalloc( sizeof(int)*LengthG) ;
    int *PKStart=(int*)mxMalloc( sizeof(int)*LengthG) ;
    int *PKEnd=(int*)mxMalloc( sizeof(int)*LengthG) ;
    int *IndexSP=(int*)mxMalloc( sizeof(int)*mxGetM(prhs[1]));
    int *IndexSQ=(int*)mxMalloc( sizeof(int)*mxGetM(prhs[0])) ;
    int MaxPk=0;
    int MaxQk=0;

    for(int i=0;i<LengthG;i++){
        // Q
        QKStart[i]=MaxQk;
        int j=MaxQk;
//         mexPrintf("%u",j);mexPrintf("\n");
        for (j=MaxQk;j<mxGetM(prhs[0]);j++){
            if (SE1[j]!=doubleUniqueG[i]){
                break;
            }
        }

        MaxQk=j;
        j=j-1;
        QKEnd[i]=min(j,mxGetM(prhs[0])-1);
        // P

        PKStart[i]=MaxPk;
       j=MaxPk;

        for (j=MaxPk;j<mxGetM(prhs[1]);j++){
            if (SE2[j]!=doubleUniqueG[i]){
                break;
            }
        }
         MaxPk=j;
        j=j-1;
        PKEnd[i]=min(j,mxGetM(prhs[1])-1);
    }
    //Q

    for(int i=0;i<LengthG;i++){
        for(int t=QKStart[i];t<=QKEnd[i];t++){
            if (t>=mxGetM(prhs[0])){
                break;
            }
          //  mexPrintf("%u",t);mexPrintf("u");mexPrintf("%f",SE1[t+2*mxGetM(prhs[0])]);mexPrintf("f");mexPrintf("%f",SE1[t+3*mxGetM(prhs[0])]);mexPrintf("======\n");
            
            double temp=Nq*(SE1[t+2*mxGetM(prhs[0])]-1)+(SE1[t+3*mxGetM(prhs[0])]);
//              mexPrintf("%u",(int)(0.5+temp));;mexPrintf("======\n");
            IndexSQ[t]=(int)(0.5+temp);
  //           mexPrintf("%u",temp);mexPrintf(",");
            int tempF=(int)mxGetM(prhs[0]);
            tempF=tempF+t;
            tempF=(int)(SE1[tempF]);
            F[t]=val[tempF-1];
        }
    }
    //P
    for(int i=0;i<LengthG;i++){
        for(int t=PKStart[i];t<=PKEnd[i];t++){
            if (t>=mxGetM(prhs[1])){
                break;
            }
            double temp=Np*(SE2[t+1*mxGetM(prhs[1])]-1)+(SE2[t+2*mxGetM(prhs[1])]);
            IndexSP[t]=(int)(0.5+temp);
        }
    }
    
    
    
    // Output
    plhs[0] = mxCreateDoubleMatrix(1,mxGetM(prhs[1]), mxREAL); //IndexSP
    plhs[1] = mxCreateDoubleMatrix(1,mxGetM(prhs[0]), mxREAL);// IndexSQ
    plhs[2] = mxCreateDoubleMatrix(1,mxGetM(prhs[0]), mxREAL); //F
    plhs[3] = mxCreateDoubleMatrix(1,LengthG, mxREAL); //PKStart
    plhs[4] = mxCreateDoubleMatrix(1,LengthG, mxREAL); //PKEnd
    plhs[5] = mxCreateDoubleMatrix(1,LengthG, mxREAL); //QKStart
    plhs[6] = mxCreateDoubleMatrix(1,LengthG, mxREAL); //QKEnd
    double *output1,*output2,*output3,*output4,*output5,*output6,*output7 ;
    output1=(double*)mxGetData(plhs[0]);//IndexSP mxGetM(prhs[1])
    output2=(double*)mxGetData(plhs[1]);//IndexSQ mxGetM(prhs[0])
    output3=(double*)mxGetData(plhs[2]);//F mxGetM(prhs[0])
    output4=(double*)mxGetData(plhs[3]); //PKStart engthG
    output5=(double*)mxGetData(plhs[4]);//PKEnd LengthG
    output6=(double*)mxGetData(plhs[5]);//QKStart LengthG
    output7=(double*)mxGetData(plhs[6]);//QKEnd LengthG
    mexPrintf("\n");
    for (int i=0;i<LengthG;i++){
        output4[i]=(double)PKStart[i]+1;
        output5[i]=(double)PKEnd[i]+1;
        output6[i]=(double)QKStart[i]+1;
        output7[i]=(double)QKEnd[i]+1;
    }
    for (int i=0;i<mxGetM(prhs[1]);i++){
        output1[i]=(double)IndexSP[i];
    }
    for (int i=0;i<mxGetM(prhs[0]);i++){
        output2[i]=(double)IndexSQ[i];
        output3[i]=(double)F[i];
    }
    mxFree(F);
    mxFree(QKStart);
    mxFree(QKEnd);
    mxFree(PKStart);
    mxFree(PKEnd);
    mxFree(IndexSP);
    mxFree(IndexSQ);
    return;
}