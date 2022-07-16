function [IsChecked,lambda,f]=FSOSChecker(P,Index,File,Lmin)
x=min(eig(P));
if x<0
    lambda=inf;
    return
end
[f,~]=DICMS2function(File);
n=f.n;
f=f-Lmin+0.5;
%% Polynomial Case
SOS=CZ_2nBig(n);
for i=1:length(P)
    for j=1:length(P)
        t=mod(Index(i,:)+Index(j,:),2);
        try
        SOS(t)=SOS(t)+P(i,j);
        catch
            disp(1)
        end
    end
end
h=f-SOS;
[~,c]=find(h);
lambda=norm(c,1);
IsChecked=0;
if lambda<0.5
    IsChecked=1;
end
end