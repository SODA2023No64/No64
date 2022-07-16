function f=CNFCharacterFunction(phi,n)
if  iscell(phi)
    f=CZ_2n(n);
    for i=1:length(phi)
        h=CNFCharacterFunction(phi{i},n);
        f=f+h;
    end
    return
end
f=CZ_2n(n);
f(2^n)=1;
for i=1:length(phi)
    h=CZ_2n(n);
    t=abs(phi(i));
    s=sign(phi(i));
    h(2^n)=1/2;
    E=zeros(1,n)+2;
    E(t)=1;
    h(Z22Qindex(E))=s/2;
    f=f.*h;
end

end


function k=Z22Qindex(g)
% {1,2}^n \mapsto N
g=(g(1)<47)*47+g;
g=char(g);
g=g(end:-1:1);
k=bin2dec(g);
k=k+1;
end