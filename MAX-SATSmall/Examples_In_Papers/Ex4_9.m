%\phi=x_1\land x_2 \land x_3\land ...\land x_{10} \lor \left(x_1 \lor x_2 \right)
%f(x)=x1/4 + x2/4 + x3/2 + x4/2 + x5/2 + x6/2 + x7/2 + x8/2 + x9/2 + x10/2 + (x1*x2)/4 + 19/4
%p_1(x)= 0.76116768833028691787490060960408+ 0.26390002447703486687657914444571*x;
%p1(f)=0.065975*x1 + 0.065975*x2 + 0.13195*x3 + 0.13195*x4 + 0.13195*x5 + 0.13195*x6 + 0.13195*x7 + 0.13195*x8 + 0.13195*x9 + 0.13195*x10 + 0.065975*x1*x2 + 2.0147
%\|p1(f)^2-f\|_{inf}>1.18
%L-1 norm |p1(f)^2-f| > 2.25
cd ..
[f,AllFai]=DICMS2function('.\Examples_In_Papers\Ex4-9.wcnf');
m=size( AllFai,1);
Lmin=1;
f=f-0.5;
varx=0.5:1:(m-Lmin-0.5);varx=varx(:);
d=1;
P_dA=varx.^[0:d];
vary=varx.^[0.5];
cvx_begin quiet
variable c(d+1,1)
minimize(norm(vary-P_dA*c,inf))
cvx_end
P1=@(x)c(1)+c(2)*x;
norm(P1(varx)-vary,inf);
g=c(1)+c(2)*f;
h=g*g-f;
[~,L1]=find(h);
L1=norm(L1,1);
Linf=CZifft(h);
Linf=max(abs(Linf(:)));
[subs,val]=find(g);
[sortval,SortInd]=sort(abs(val),'descend');
subs=subs(SortInd);
[FP,P,ti,SOSP,lambda,cvx_cputime]=PolySOSWithL_1Error(f,subs(1:end-1));
cd Examples_In_Papers
disp(L1);
disp(Linf);
disp(lambda);
norm(P1(varx).^2-varx,inf)