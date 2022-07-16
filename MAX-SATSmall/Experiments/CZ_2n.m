classdef CZ_2n %< handle
    % 结构：CZ_2n.n 就表示Z_2^n的n，
    % fhat=CZ_2n.c fhat(i)就是 CZ_2n的系数 既： CZ_2n.c(i)就是Qindex2Z2(i,CZ_2n.n)项的系数
    properties
        n
        c
    end
    methods
        %% 构造函数
        function obj = CZ_2n(n)
            if isnumeric(n)&&(length(n)==1)% 是数字
                obj.n=n;
                obj.c= containers.Map('KeyType',  'double', 'ValueType', 'any');%空函数
                return
            end
            %             else  % 是映射
            %                 obj.n=n;
            %                 obj.c=c;
            if isnumeric(n)&&(isvector(n)==1)&&(length(n)>1)%是sparse向量或者full向量
                m=log2(length(n));
                m=ceil(m);
                Index=find(n);
                c=containers.Map(Index,full(n(Index)));
                obj.n=m;
                obj.c=c;
                return
            end
        end
        %% Get/Set
        function  x = get(f,i)
            %得到 函数的系数
            if length(i)>1
                i=Z22Qindex(i);
            end
            if isKey(f.c, i)
                x=f.c(i);
            else
                x=0;
            end
        end
        
        function  f = set(f,i,x)
            f.c(i)=x;
        end
        %% 赋值
        function x = subsref(this,s)
            
            if strcmp(s.subs,'c')
                x=this.c;
                return
            end
            if strcmp(s.subs,'n')
                x=this.n;
                return
            end
            C=this.c;
            s=s.subs;
            s=s{1};
            if isrow(s)&&~isscalar(s)
                s=mod(s,2);
                s(s==0)=2;
                s=Z22Qindex(s);
            end
            if isKey(this.c,s)
                x=C(s);
            else
                x=0;
            end
        end
        
        
        function this = subsasgn(this,s,b)
            s=s.subs;
            s=s{1};
            if isrow(s)&&~isscalar(s)
                s=mod(s,2);
                s(s==0)=2;
                s=Z22Qindex(s);
            end
            t= this.c;
            t(s)=b;
            this.c=t;
        end
        %% 加法
        function h=plus(f,g)
            if isnumeric(f)&&length(f)==1
                h=g;
                set(h,2^(g.n),get(h,2^(g.n))+f);
                return 
            end
            if isnumeric(g)&&length(g)==1
                h=plus(g,f);
                return
            end
            if f.n~=g.n
                error('the addition are not on same group')
            end
            h=CZ_2n(f.n);
            A=f.c;
            B=g.c;
            Index=[cell2mat(A.keys()),cell2mat(B.keys())];
            Index=unique(Index);
            for i=(Index)
                if get(f,i)+get(g,i)~=0
                    set(h,i,get(f,i)+get(g,i));
                end
            end
        end
        
        function L=length(this)
            L=length(this.c);
        end
        
        %% 减法
        function h = minus(f,g)
            if isnumeric(f)&&length(f)==1
                h=minus(CZ_2n(g.n),g)+f;
                return
            end
            if isnumeric(g)&&length(g)==1
                h=plus(f,-g);
                return
            end
            if f.n~=g.n
                error('the addition are not on same group')
            end
            h=CZ_2n(f.n);
            A=f.c;
            B=g.c;
            Index=[cell2mat(A.keys()),cell2mat(B.keys())];
            Index=unique(Index);
            for i=(Index)
                if get(f,i)-get(g,i)~=0
                    set(h,i,get(f,i)-get(g,i));
                end
            end
        end
        %% 乘法
        function h = mtimes(f,g)
            %C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
            if isnumeric(f)&&length(f)==1% is  scalar
                h=CZ_2n(g.n);
                for i=cell2mat(keys(g.c))
                    set(h,i,f*get(g,i));
                end
                return
            end
            if isnumeric(g)&&length(g)==1% is  scalar
                h=g*f;
                return
            end
            N=f.n;
            h=CZ_2n(N);
            A=f.c;
            B=g.c;
            Index1=cell2mat(A.keys());
            Index2=cell2mat(B.keys());
            % QindexAdd(a,b,n)=2^n-bitxor(x-1,y-1)
            %T=2^(f.n)-bitxor(Index1(:)-1,Index2(:)'-1);
            for i=Index1
                
                for j=Index2
                    x=Qindex2Z2(i,N);
                    y=Qindex2Z2(j,N);
                    t=mod(x+y,2);
                    t(t==0)=2;
                    t=Z22Qindex(t);
                    %disp([t,x,y])
                    
                    set(h,t,get(h,t)+get(f,x)*get(g,y));
                end
            end
        end
        %% 乘法2 .*
        function h=times(f,g)
            %C = TIMES(A,B) is called for the syntax 'A .* B'
            if isnumeric(f)&&length(f)==1% is  scalar
                h=CZ_2n(g.n);
                for i=cell2mat(keys(g.c))
                    set(h,i,f*get(g,i));
                end
                return
            end
            if isnumeric(g)&&length(g)==1% is  scalar
                h=g*f;
                return
            end
            N=f.n;
            h=CZ_2n(N);
            [Indexf,fhat]=find(f);
            [Indexg,ghat]=find(g);
            H=fhat(:)*ghat(:).';
            Bp=bitxor(2^N-Indexf(:),2^N-Indexg(:)');
            Bp=2^N-Bp;
            BpIndex=unique(Bp(:));
            % %             for i=1:length(BpIndex)
            % %                 if sum(sum((Bp==BpIndex(i)).*H))~=0
            % %                    % disp(i./length(BpIndex))
            % %                     set(h,BpIndex(i),sum(sum((Bp==BpIndex(i)).*H)));
            % %                 end
            % %             end
            for i=1:size(Bp,1)
                for j=1:size(Bp,2)
                    set(h,Bp(i,j),get(h,Bp(i,j))+H(i,j));
                end
            end
        end
        
        
        
        
        %% Display
        function disp(f)
            S=sym(f);
            disp(S)
        end
        %% sym多项式化
        function S=sym(f)
            x=sym('x',[f.n,1]);
            S=sym(0);
            A=f.c;
            Index=cell2mat(A.keys());
            for i=Index
                t=Qindex2Z2(i,f.n);
                t=mod(t,2);
                S=S+get(f,i)*prod(x(t>0));
            end
        end
        
        %% Find
        function [subs,vals] = find(f)
            A=f.c;
            subs=keys(A);
            subs=cell2mat(subs);
            vals=values(A);
            vals=cell2mat(vals);
        end
        
        %% polyval
        function h=polyval(P,this)
            %combine P(f)
            % where P is a vector as P(x)=P(1)*x^(N-1)+P(2)*x^(N-2)...+P(N)
            h=CZ_2n(this.n);
            g=this;
            P=P(end:-1:1);
            set(h,2^(this.n),P(1));
            for i=2:length(P)
                h=h+P(i)*g;
                if (i<length(P))
                    g=g.*this;
                end
            end
        end
        %% CZ_2n to Sparse
        function s=sparse(f)
            [subs,vals] = find(f);
            s=sparse(subs',ones(length(f),1),vals,2^(f.n),1);
        end
        %% IFFT
        function V=CZifft(f)
            s=sparse(f);
            s=full(s);
            V=ifftn(reshape(flipud(s),2*ones(1,f.n)))*2^(f.n);
            V=reshape(flipud(V(:)),2*ones(1,f.n));
        end
        %         function f=simplify(f)
        %             %remove values equals to 0
        %
        %         end
    end
       
    
end



function I=Qindex2Z2(i,siz)
%Qindex to Z_2^n : N \mapsto n*1 vectors {1,2}^n
% s.t.  T is a 2*2*2...*2 double array, T(ind2Z2n(siz,ndx))=T(ndx)

if length(siz)>1
    siz=length(siz);
end
T=dec2bin(i-1,siz);
I=(T(end:-1:1)-47);
end


function k=Z22Qindex(g)
% {1,2}^n \mapsto N
g=(g(1)<47)*47+g;
g=char(g);
g=g(end:-1:1);
k=bin2dec(g);
k=k+1;
end