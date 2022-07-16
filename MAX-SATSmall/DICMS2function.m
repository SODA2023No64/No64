function [f,AllFai]=DICMS2function(FileName)
fid=fopen(FileName,'r');
data=[];
tline=fgetl(fid);
NM_data=tline;pat='\d+';
NM_data=regexp(NM_data,pat,'match');
n=NM_data{1};n=str2num(n);
while 1
    tline=fgetl(fid);
    if~ischar(tline)
        break;
    end
    tline=str2num(tline);
    ld=size(data,1);
    for k=1:length(tline)
        data(ld+1,k)=tline(k);
    end
end
fclose(fid);
data=data(:,2:end-1);
AllFai=zeros(size(data,1),n);
for i=1:size(data,1)
    t=data(i,:);
    t=t(t~=0);
    AllFai(i,abs(t))=sign(t);
end
D=AllFai2CNF(AllFai);
f=CNFCharacterFunction(D,n);
end