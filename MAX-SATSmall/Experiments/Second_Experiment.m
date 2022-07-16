%The second experiment
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
    IndexQ=subs;
    Identy=2^f.n;
    for j=1:4
        [FP{i,j},P{i,j},ti{i,j},SOSPe{i,j},lambda{i,j},cvx_cputime{i,j}]=PolySOSWithL_1Error(f,IndexQ(1:end*(j/(j+1))),1);
    end
    j=5;
    [FP{i,j},P{i,j},ti{i,j},SOSPe{i,j},lambda{i,j},cvx_cputime{i,j}]=PolySOSWithL_1Error(f,IndexQ(1:end),1);
    clear FP P ti SOSPe lambda cvx_cputime4;
end