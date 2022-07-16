function [FP,P,ti,SOSP,lambda,cvx_cputime]= PolySOSWithL_1Error(f,IndexP,ComputeFlag)
%IndexP£º Support of polynomial
% ComputeFlag: if exists then only compute lambda and P and Q
cvx_begin;cvx_end;cvx_clear;
IndexP=IndexP(:);
Indexf=find(f);Indexf=Indexf(:);
disp('Start Compute')
n=f.n;
ti=cputime;
F=NaN;Q=NaN;
cvx_status='0';
cvx_precision default
%%  QindexAdd(a,b,n)=2^n-bitxor(x-1,y-1)
%% condition : sum( f(a)*Q(b,c),abc=t)= sum( P(d,e), de=t) for all t.
Identy=2^(f.n);
IndexQ=Identy;
[SparseEq1,SparseEq2]=RationalCondition(Identy,IndexP,Indexf,n);
SE1=sortrows(SparseEq1,1);
SE2=sortrows(SparseEq2,1);
UniqueG=unique([SE1(:,1);SE2(:,1)]);UniqueG=UniqueG(:);UniqueG=sort(UniqueG);
GNq=1;
GNp=length(IndexP);
[subs,val]=find(f);
[IndexSP,IndexSQ,F,PKStart,PKEnd,QKStart,QKEnd]=AfterCondition(SE1,SE2,UniqueG,subs,val,GNp,GNq);
clear SE1 SE2
Aq=sparse(length(UniqueG(:)),GNq*GNq);
Ap=sparse(length(UniqueG(:)),GNp*GNp);
for i=1:length(UniqueG(:))
    Aq(i,IndexSQ(QKStart(i):QKEnd(i)))=F(QKStart(i):QKEnd(i));
    Ap(i,IndexSP(PKStart(i):PKEnd(i)))=1;
end
try
    cvx_begin sdp
    variable P(GNp,GNp)
    minimize(norm(Aq-Ap*P(:),1))
    P >= 0;
    cvx_condition_end=cputime;
    cvx_end
    lambda=cvx_optval;
    Q=1;
catch
    if strcmp(cvx_status,'Error')
        SOSP=inf;
        FP=inf;
        lambda=inf;
        cvx_cputime=inf;
        ti=inf;
        return
    end
end
ti=cputime-ti;
cvx_time_end=cputime;
disp(['sparsity:' num2str(length(find(diag(Q)))),',',num2str(length(find(diag(P))))])
cvx_cputime=cvx_time_end-cvx_condition_end;
if ((nargin==3)&&(ComputeFlag==1))
    SOSP=inf;
    FP=inf;
    return
end
[LQ,DQ]=ldl(full(Q));
LQ=sqrt(DQ)*LQ';
FQ=cell(length(IndexQ),1);
for i=1:size(LQ,1)
    FQ{i}=CZ_2n(n);
    for j=1:length(IndexQ)
        FQ{i}(IndexQ(j))=LQ(i,j);
    end
end

[LP,DP]=ldl(full(P));
LP=sqrt(DP)*LP';
FP=cell(length(IndexP),1);
for i=1:size(LP,1)
    FP{i}=CZ_2n(n);
    for j=1:length(IndexP)
        FP{i}(IndexP(j))=real(LP(i,j));
    end
end
SOSP=CZ_2n(n);
for i=1:length(FP)
    SOSP=SOSP+FP{i}.*FP{i};
end

end
