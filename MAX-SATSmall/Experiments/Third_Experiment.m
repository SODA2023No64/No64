%The third experiment
clear
load('RandomResult.mat')
cd ..
for i=1:70
    FileName=string('.\output\')+File{i};
    cd Experiments
    [f,AllFai]=DICMS2function(FileName);
    cd ..
    f=f-MS(i)+0.5;
    [subs,val]=find(f);
    [sortval,SortInd]=sort(abs(val),'descend');
    subs=subs(SortInd);
    S=length(subs);
    S=floor(S/2);
    Identy=2^f.n;
    f1=f;
    f2=f.*f;
    for d=1:2
        mind=sqrt(length(subs));
        mind=floor(mind/2);%minimum unfeasible $s$
        mind=max(mind,1);%maximum feasible $s$
        maxd=inf;
        m=size(AllFai,1);
        varx=0.5:1:(m-MS(i)-0.5);varx=varx(:);
        P_dA=varx.^[0:d];
        vary=varx.^[0.5];
        cvx_begin quiet
        variable c(d+1,1)
        minimize(norm(vary-P_dA*c,inf))
        cvx_end
        c(5)=0;
        g=c(1)+c(2)*f1+c(3)*f2;
        [subs,val]=find(g);
        [sortval,SortInd]=sort(abs(val),'descend');
        subs=subs(SortInd);
        %% Enlarge stage:
        lambda=inf;
        
        while (lambda>=0.5)
            mind=min(2*mind,length(subs));
            [~,~,~,~,lambda,~]=PolySOSWithL_1Error(f,subs(1:min(2*mind,length(subs))),1);
        end
        maxd=min(2*mind,length(subs));
        %% 1/2 stage:
        while (maxd>mind+1)
            s=floor((maxd+mind)*0.5);
            [~,~,~,~,lambda,~]=PolySOSWithL_1Error(f,subs(1:s),1);
            pause(3);
            if lambda>=0.5
                mind=s;
            else
                maxd=s;
            end
        end
        Sparse(i,d)=maxd;
    end
    clear FP P ti SOSPe cvx_cputime4 erre err;

end
