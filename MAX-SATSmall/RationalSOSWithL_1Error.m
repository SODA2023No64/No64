function [FP,FQ,P,Q,ti,SOSP,SOSQ,lambda,cvx_cputime]= RationalSOSWithL_1Error(f,IndexP,IndexQ,ComputeFlag)
%IndexP£ºsupp(h_j) £¬ IndexQ£ºsupp(g_i)
% ComputeFlag: if exists then only compute lambda and P and Q
cvx_begin;cvx_end;cvx_clear;
IndexP=IndexP(:);
IndexQ=IndexQ(:);
Indexf=find(f);Indexf=Indexf(:);
disp('Start Compute')
n=f.n;
ti=cputime;
F=NaN;Q=NaN;
cvx_status='0';
cvx_precision default
%%  QindexAdd(a,b,n)=2^n-bitxor(x-1,y-1)
%% condition : sum( f(a)*Q(b,c),abc=t)= sum( P(d,e), de=t) for all t.
[SparseEq1,SparseEq2]=RationalCondition(IndexQ,IndexP,Indexf,n);
SE1=sortrows(SparseEq1,1);
SE2=sortrows(SparseEq2,1);
UniqueG=unique([SE1(:,1);SE2(:,1)]);UniqueG=UniqueG(:);UniqueG=sort(UniqueG);
GNq=length(IndexQ);
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
    variable Q(GNq,GNq)
    minimize(norm(Aq*Q(:)-Ap*P(:),1))
    cvx_time_start=cputime;
    Q >= 1/GNq*eye(GNq)
    P >= 0;
    cvx_condition_end=cputime;
    cvx_end
    lambda=cvx_optval;
catch
    if strcmp(cvx_status,'Error')
        SOSP=inf;
        SOSQ=inf;
        FP=inf;
        FQ=inf;
        lambda=inf;
        cvx_cputime=inf;
        ti=inf;
        return
    end
end
ti=cputime-ti;
cvx_time_end=cputime;
%cvx_solver_settings('MSK_IPAR_INTPNT_MULTI_THREAD',0);
disp(['sparsity:' num2str(length(find(diag(Q)))),',',num2str(length(find(diag(P))))])
cvx_cputime=cvx_time_end-cvx_condition_end;
if (nargin==3)||(ComputeFlag==0)
    SOSP=inf;
    SOSQ=inf;
    FP=inf;
    FQ=inf;
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
SOSQ=CZ_2n(n);
for i=1:length(FQ)
    SOSQ=SOSQ+FQ{i}.*FQ{i};
end

end
