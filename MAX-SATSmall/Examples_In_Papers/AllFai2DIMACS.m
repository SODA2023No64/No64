function D=AllFai2DIMACS(AllFai,FileName,txtflag)
%input: AllFai
%output:DIMACS format 
%if input FileName exists,then save it as FileName.txt

    n=size(AllFai,2);% number of variables
    m=size(AllFai,1);% number of clauses
    D{1}=['p wcnf ' num2str(n) ,' ',num2str(m) ,' ',num2str(m+1) ];
    for i=2:m+1
        x=AllFai(i-1,:);
        t=[1,-find(x==-1),find(x==1),0];
        D{i}=num2str(t);
    end
if nargin==2
fid=fopen([FileName,'.txt'],'w');
for i=1:length(D)
fprintf(fid,D{i});
fprintf(fid,'\n');
end
fclose(fid);
end
if nargin==3
fid=fopen([FileName],'w');
for i=1:length(D)
fprintf(fid,D{i});
fprintf(fid,'\n');
end
fclose(fid);
end
end