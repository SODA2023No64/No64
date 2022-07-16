function f=CNFCharacterFunction2(phi,n,w)
if  iscell(phi)
    if nargin==2
        w=ones(length(phi),1);
    end
end
f=CZ_2nBig(n);
K=0;
m=0;
for i=1:length(phi)
    K=max(K,length(phi{i}));
    m=m+2^length(phi{i});
end
NcK=dec2bin(0:(2^K-1),K);
NcK=NcK-48;
A=-1+zeros(m,n);
b=zeros(m,1);
s=0;
for i=1:length(phi)
    d=length(phi{i});
    t=abs(phi{i});
    sig=(phi{i})<0;
    sig=double(sig)';
    A(s+1:s+2^d,:)=0;
    A(s+1:s+2^d,t)=NcK(1:2^d,(K-d+1):K);
    bt=NcK(1:2^d,(K-d+1):K)*sig;
    b(s+1:s+2^d)=w(i)*(-1).^bt/2^d;
    %1/(2^d)
    s=s+2^d;
end
[SA,IA]=sortrows(A);%SA=A(IA,:)
Sb=b(IA);
t=SA(1,:);
sumb=0;
for i=1:size(SA,1)
    if all(SA(i,:)==t)
        sumb=sumb+Sb(i);
    else
        if sumb~=0
            f(t)=sumb;
        end
        t=SA(i,:);
        sumb=Sb(i);
    end
end
        if sumb~=0
            f(t)=sumb;
        end
end