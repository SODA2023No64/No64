%The first expperiment
% This is also an exapmle to show how to use FSOSBuilder
% Lmin is th minimal value of \phi
clear
load('RandomResult.mat')
Lmin=MS;
cd ..
for i=1:70
FileName=string('.\Experiments\output\')+File{i};
T_Start(i)=cputime;
[P,lambda(i),Index,f,F,SOS,err]= FSOSBuilder(FileName,Lmin(i),0,0);
%[P,lambda(i),Index,f,F,SOS,err]= FSOSBuilder(FileName,MS(i),1,0); For Rational FSOS
T_End(i)=cputime;
T(i)=T_End(i)-T_Start(i);
end
cd Experiments