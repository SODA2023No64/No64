function [lambda,f,p]=FSOSChecker(P,Index,File,Lmin)
x=min(eig(P));
if x<0
    lambda=inf;
    return
end
[f,~]=DICMS2function(File);
n=f.n;
f=f-Lmin+0.5;
if (~iscell(Index))
    %% Polynomial Case
    ProdIndex=2^n-bitxor(2^n-Index(:),2^n-Index(:)');
    p=CZ_2n(n);
    for i=1:length(ProdIndex)
        for j=1:length(ProdIndex)
            p(ProdIndex(i,j))=p(ProdIndex(i,j))+P(i,j);
        end
    end
    h=f-p;
    [~,c]=find(h);
    lambda=norm(c,1);
else
    %% Rational Case
    ProdIndexP=bitxor(2^n-Index{1}(:),2^n-Index{1}(:)');
    ProdIndexQ=bitxor(2^n-Index{2}(:),2^n-Index{2}(:)');
    p=CZ_2n(n);q=CZ_2n(n);
    for i=1:length(ProdIndexP)
        for j=1:length(ProdIndexP)
            p(ProdIndexP(i,j))=p(ProdIndexP(i,j))+P{1}(i,j);
        end
    end
    for i=1:length(ProdIndexQ)
        for j=1:length(ProdIndexQ)
            q(ProdIndexQ(i,j))=q(ProdIndexQ(i,j))+P{2}(i,j);
        end
    end
    h=p-q*f;
    [~,c]=find(h);
    lambda=norm(c,1);
end

end