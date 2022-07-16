function D=AllFai2CNF(AllFai)
    n=size(AllFai,2);% number of variables
    m=size(AllFai,1);% number of clauses
    for i=1:m
        x=AllFai(i,:);
        t=[-find(x==-1),find(x==1)];
        D{i}=t;
    end
end