function [P,ti,lambda,cvx_cputime]= PolySOSWithL_1Error(f,IndexP)
%IndexP£º Support of polynomial In form of  m*n 0-1 matrix s.t. IndexP(i,:)
%is an element of IndexP
cvx_begin;cvx_end;cvx_clear;
Indexf=find(f);
disp('Start Compute')
n=f.n;
ti=cputime;
% F=NaN;Q=NaN;
cvx_status='0';
cvx_precision default

Identy=zeros(1,n);
[SparseEq1,SparseEq2]=RationalCondition(Identy,IndexP,Indexf,n);
SEQ=sortrows(SparseEq1,1:n);%SE1
SEP=sortrows(SparseEq2,1:n);%SE2
UniqueG=unique([SEQ(:,1:n);SEP(:,1:n)],'rows');
GNq=1;
GNp=size(IndexP,1);
[subs,val]=find(f);
[IndexSP,IndexSQ,F,PKStart,PKEnd,QKStart,QKEnd]=AfterCondition(SEQ,SEP,UniqueG,subs,val,GNp,GNq);
%clear SEP SEQ
Aq=sparse(size(UniqueG,1),GNq*GNq);
Ap=sparse(size(UniqueG,1),GNp*GNp);
for i=1:(size(UniqueG,1))
    Aq(i,IndexSQ(QKStart(i):QKEnd(i)))=F(QKStart(i):QKEnd(i));
    Ap(i,IndexSP(PKStart(i):PKEnd(i)))=1;
end
% try

%Method 1
cvx_begin sdp
variable P(GNp,GNp) symmetric
%minimize(norm(Aq-Ap*P(:),1))
norm(Aq-Ap*P(:),1)<=0.5;
P >= 0;
cvx_condition_end=cputime;
cvx_end

% % % %Method 2
% cvx_begin sdp
% variable P(GNp,GNp)
% %variable Lam(size(UniqueG,1),1)
% B=zeros(size(UniqueG,1),1);
% B=cvx(B);
% for i=1:size(UniqueG,1)
%     try
%     B(i,1)=Aq(i)-sum([0;P(IndexSP(PKStart(i):PKEnd(i)))]);
%     catch 
%         disp(1)
%     end
% end
% norm(B(:),1)<=0.5;
% P >= 0;
% cvx_condition_end=cputime;
% cvx_end

% Method 3
% % Aq(2:end,:)=Aq(2:end,:)/2;
% % Ap=sparse(size(UniqueG,1),GNp*GNp);
% % for i=1:(size(UniqueG,1))
% %     T=SEP(PKStart(i):PKEnd(i),n+1:end);
% %     T=T( T(:,1)<T(:,2),:);
% %     T=T*[GNp;1]-GNp;
% %     Ap(i,T)=1;
% % end
% % cvx_begin sdp
% % variable P(GNp,GNp)
% % norm(Aq-Ap*P(:),1)<=0.5;
% % P >= 0;
% % cvx_condition_end=cputime;
% % cvx_end

lambda=cvx_optval;
lambda=norm(Aq-Ap*P(:),1);
Q=1;
% catch
if strcmp(cvx_status,'Error')
    lambda=inf;
    cvx_cputime=inf;
    ti=inf;
    return
end
% end
ti=cputime-ti;
cvx_time_end=cputime;
disp(['sparsity:' num2str(length(find(diag(Q)))),',',num2str(length(find(diag(P))))])
cvx_cputime=cvx_time_end-cvx_condition_end;
end


function [SparseEq1,SparseEq2]=RationalCondition(IndexQ,IndexP,Indexf,n)
n=size(IndexP,2);
Lp=size(IndexP,1);
Lq=size(IndexQ,1);
SparseEq2=zeros(Lp*Lp,n+2);
for i=1:Lp
    for j=1:Lp
        tempi=Lp*(i-1)+j;
        t=mod(IndexP(i,:)+IndexP(j,:),2);
        SparseEq2(tempi,:)=[t,i,j];
    end
end
SparseEq1=zeros(size(Indexf,1)*Lq*Lq,n+3);
for i=1:size(Indexf,1)
    for j=1:Lq
        for k=j:Lq
            tempi=Lq*(j-1)+k+Lq*Lq*(i-1);
%             disp([k,j,i,tempi])
            t=mod(IndexQ(j,:)+IndexQ(k,:)+Indexf(i,:),2);
            SparseEq1(tempi,:)=[t,i,j,k];
        end
    end
end
end

%@(x,y)(n-1)*x+y-(n-1)
% for i=1:Lp
%     for j=i:Lp
%         tempi=(Lp-1)*i+j-(Lp-1);
%         t=mod(IndexP(i,:)+IndexP(j,:),2);
%         SparseEq2(tempi,:)=[t,i,j];
%     end
% end
% SparseEq1=zeros(size(Indexf,1)*(Lq+1)*Lq/2,n+3);
% for i=1:size(Indexf,1)
%     for j=1:Lq
%         for k=j:Lq
%             tempi=(Lp-1)*j+k-(Lp-1)+(i-1)*Lq*(Lq+1)/2;
%             disp([k,j,i,tempi])
%             t=mod(IndexQ(j,:)+IndexQ(k,:)+Indexf(i,:),2);
%             SparseEq1(tempi,:)=[t,i,j,k];
%         end
%     end
% end



function [IndexSP,IndexSQ,F,PKStart,PKEnd,QKStart,QKEnd]=AfterCondition(SEQ,SEP,UniqueG,subs,val,GNp,GNq)
G=size(UniqueG,1);
PKStart=zeros(G,1);
PKEnd=zeros(G,1);
QKStart=zeros(G,1);
QKEnd=zeros(G,1);
Np=size(SEP,1);
Nq=size(SEQ,1);
F=zeros(Nq,1);
n=size(UniqueG,2);
% sum( P(IndexP(PKStart(i):PKEnd(i))))==sum( Q(IndexQ(QKStart(i):QKEnd(i))).*F(QKStart(i):QKEnd(i)))
MaxPk=1;
MaxQk=1;
for i=1:size(UniqueG,1)
    QKStart(i)=MaxQk;
    j=MaxQk;
    LastEqualFlag=0;
    for j=MaxQk:size(SEQ,1)
        if ~all(SEQ(j,1:n)==UniqueG(i,:))
            LastEqualFlag=1;
            break;
        end
    end
    MaxQk=j;
    j=j-LastEqualFlag;
    QKEnd(i)=min(j,size(SEQ,1));
    PKStart(i,:)=MaxPk;
    LastEqualFlag=0;
    for k=MaxPk:size(SEP,1)%(j=MaxPk;j<mxGetM(prhs[1]);j++){
        if  ~all(SEP(k,1:n)==UniqueG(i,:))%(SE2[j]!=doubleUniqueG[i]){
            LastEqualFlag=1;
            break;
        end
    end
    MaxPk=k;
    k=k-LastEqualFlag;
    PKEnd(i)=min(k,Np);
    IndexSP2=SEP(:,n+1:end);
    IndexSP=IndexSP2(:,1)+(IndexSP2(:,2)-1)*GNp;
    IndexSQ2=SEQ(:,n+2:end);
    IndexSQ=IndexSQ2(:,1)+(IndexSQ2(:,2)-1)*GNq;
end
F=zeros(length(IndexSQ),1);
for i=1:G
    F(QKStart(i):QKEnd(i))=val(SEQ(QKStart(i):QKEnd(i),n+1));
end
end