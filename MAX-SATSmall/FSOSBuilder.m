function [P,lambda,Index,f,F,SOS,err]= FSOSBuilder(File,Lmin,Rational_Flag,Check_Flag)
%% This is the main function of FSOS verificiation
% input:
%  File: file name of wcnf phi, file should write in wcnf form
%  Lmin: min_I cost_I (\phi)
%  Rational_Flag=1 then compute the Rational FSOS, otherwise Polynomial (default 0)
%  iFFT_Check_Flag=1 then compute pointwise error of FSOS, otherwise dont compute it(default 0)
% [P,lambda,Index,f]= FSOSBuilder(File,Lmin,0,0) computes the polynomial
% FSOS of f:=f_{\phi}-Lmin+1/2 of given CNF in File, P is the Gram matrix, lambda is the Fourier L_1 error of SDP problem.
% Index is the label of P.
%
% [P,lambda,Index,f]= FSOSBuilder(File,Lmin,1,0)  computes the rational
% FSOS of f:=f_{\phi}-Lmin+1/2 of given CNF in File, P{1} and P{2} are the Gram matrix of numerator and denominator,
% lambda is the Fourier L_1 error of SDP problem. Index{1} and Index{2} are
% the label of  Gram matrix of numerator and denominator.
%
% [P,lambda,Index,f,F,SOS,err]= FSOSBuilder(File,Lmin,0,1) and [F,P,ti,SOSP,lambda,err]= FSOSBuilder(File,Lmin,1,1)
% not only computes the polynomial/rational Gram matrix of FSOS but  also
% compute the form of FSOS. here F is a cell with  F{:,1} is the terms of
% numerator (and F{:,2} is the terms of enominator if compute rational FSOS).
% SOS{1} and SOS{2} are the expansion of numerator and denominator
% (Polynomial case only exists SOS{1}) and  also use inverse discrete Fourier transform to compute the maximal pointwise error of FSOS and f_{\phi}-Lmin+1/2.

%% Start Computation
err=nan;
[f,AllFai]=DICMS2function(File);
n=f.n;
m=size(AllFai,1);
f=f-Lmin+0.5;
if nargin==2
    Rational_Flag=0;Check_Flag=0;
end
if nargin==3
    Check_Flag=0;
end
d=0;
lambda=10;
fd{1}=CZ_2n(n)+1;
while lambda>=0.5
    d=d+1;
    %%     compute P_d and P_d(f)
    varx=0.5:1:(m-Lmin-0.5);varx=varx(:);
    P_dA=varx.^[0:d];
    vary=varx.^[0.5];
    cvx_begin quiet
    variable c(d+1,1)
    minimize(norm(vary-P_dA*c,inf))
    cvx_end
    Approach_d=cvx_optval;
    P_d_f=CZ_2n(n);
    fd{d+1}=fd{d}*f;
    for i=1:(d+1)
        P_d_f=P_d_f+c(i)*fd{i};
    end
    if Rational_Flag~=1
        %% compute polynomial FSOS
        for IndexRate=[1/3,1/2,2/3,3/4,4/5,1] %\pho
            [subs,val]=find(P_d_f);
            [sortval,SortInd]=sort(abs(val),'descend');
            subs=subs(SortInd);
            Index=subs(1:end*IndexRate);
            [F,P,ti,SOS,lambda,cvx_cputime]=PolySOSWithL_1Error(f,Index,1-Check_Flag);
            if lambda<0.5
                if Check_Flag==1
                    erre=CZifft(SOS)-CZifft(f);
                    err=max(abs(erre(:)));
                end
                break;
            end
        end
    else
        %% rational FSOS
        for IndexRate=[1/3,1/2,2/3,3/4,4/5,1] %\pho
            [subs,val]=find(P_d_f);
            [sortval,SortInd]=sort(abs(val),'descend');
            subs=subs(SortInd);
            IndexQ=subs(1:end*IndexRate);
            IndexP=2^(f.n)-bitxor(subs(:)'-1,IndexQ(:)-1);
            IndexP=unique(IndexP(:));
            [FP,FQ,P,Q,ti,SOSP,SOSQ,lambda]= RationalSOSWithL_1Error(f,IndexP,IndexQ,Check_Flag);
            if lambda<0.5
                if Check_Flag==1
                    erre=CZifft(SOSP)-CZifft(f).*CZifft(Q);
                    err=max(abs(erre(:)));
                    Index{1}=IndexP;Index{2}=IndexQ;
                    Pt=P;
                    clear P
                    P{1}=Pt;P{2}=Q;
                    SOS{1}=SOSP;
                    SOS{2}=SOSQ;
                    F{1}=FP;
                    F{2}=FQ;
                else
                    F=inf;
                    err=inf;
                    SOS=inf;
                    Index{1}=IndexP;
                    Index{2}=IndexQ;
                    Pt=P;
                    clear P
                    P{1}=Pt;P{2}=Q;
                end
                break;
            end
        end
    end
    
end

end