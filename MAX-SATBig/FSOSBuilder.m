function [P,lambda,Index,f]= FSOSBuilder(File,Lmin)
%% This is the main function of FSOS  verificiation (ONLY POLORMIAL!)
% input:
%  File: file name of wcnf phi, file should write in wcnf form
%  Lmin: min_I cost_I (\phi)

% [P,lambda,Index,f]= FSOSBuilder(File,Lmin)  computes the polynomial
% FSOS of f:=f_{\phi}-Lmin+1/2 of given CNF in File, P is the Gram
% matrix, lambda is the opt value of SDP porblem and Index is the label of
% matrix P

[f,w]=DICMS2function(File);
m=sum(w);
f=f-(Lmin-0.5);
lambda=inf;
d=0;
n=f.n;
fd{1}=CZ_2nBig(n)+1;
while (lambda >=0.5)
    d=d+1;
    varx=0.5:1:(m-Lmin-0.5);varx=varx(:);
    P_dA=varx.^[0:d];
    vary=varx.^[0.5];
    cvx_begin quiet
    variable c(d+1,1)
    minimize(norm(vary-P_dA*c,inf))
    cvx_end
    Approach_d=cvx_optval;
    P_d_f=CZ_2nBig(n);
    fd{d+1}=fd{d}*f;
    for i=1:(d+1)
        P_d_f=P_d_f+c(i)*fd{i};
    end
    %% compute polynomial FSOS
    [subs,val]=find(P_d_f);
    [~,SortInd]=sort(abs(val),'descend');
    subs=subs(SortInd,:);
    for r=[0.5,0.8,1]
        Index=subs(1:round(r*length(val)),:);
        [P,~,lambda,~]=PolySOSWithL_1Error(f,Index);
        if lambda<0.5
            break;
        end
    end
end
end