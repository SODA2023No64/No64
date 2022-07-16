% The Fourth Experiment 
Index=1:10;
 Lmin=[14 12 13 18 13 24 28 25 25 29] ;
X=dir('.\BigExams');
k=0;
disp(Index)
for i=Index
    disp(['Compute: ', num2str(i)])
    N=['.\BigExams\',num2str(i) '.wcnf'];
    k=k+1;
    File{k+1}=N;
    L=Lmin(i);
    T(i)=cputime;
    [P,lambda,Index,f]= FSOSBuilder(N,L);
    T(i)=cputime-T(i);
end