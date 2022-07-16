function [f,w]=DICMS2function2(FileName)
fid=fopen(FileName,'r');
data=[];
tline=fgetl(fid);
while( ~strcmp(tline(1:min(6,length(tline))),'p wcnf'))
    tline=fgetl(fid);
end
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
w=data(:,1);
data=data(:,2:end-1);
for i=1:length(w)
    ld=data(i,:);
    ld=ld(ld~=0);
    D{i}=ld;
end
f=CNFCharacterFunction2(D,n,w);
end