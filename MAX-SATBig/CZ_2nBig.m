classdef CZ_2nBig %< handle
    % 结构：CZ_2n.n 就表示Z_2^n的n，
    % fhat=CZ_2n.c fhat(i)就是 CZ_2n的系数 既： CZ_2n.c(i)就是Qindex2Z2(i,CZ_2n.n)项的系数
    properties
        n
        c
    end
    methods
        %% 构造函数
        function obj = CZ_2nBig(n)
            if isnumeric(n)&&(length(n)==1)% 是数字
                obj.n=n;
                obj.c=containers.Map('KeyType',  'char', 'ValueType', 'any');
                return
            end
        end
        %% Get/Set
        function  x = get(f,t)
            %得到 函数的系数
            if ~ischar(t)
                t=char(t+48);
            end
            if length(t(:))==f.n
                t=t(:)';
                if isKey(f.c,t)
                    x=f.c(t);
                else
                    x=0;
                end
            end
        end
        
        function  f = set(f,t,x)
            if ~ischar(t)
                t=mod(t,2);
                t=char(t+48);
            end
            f.c(t)=x;
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
            s=s.subs;
            t=s{1};
            if ~ischar(t)
                t=mod(t,2);
                t=char(t+48);
            end
            if isKey(this.c,t)
                x=get(this,t);
            else
                x=0;
            end
        end
        
        function this = subsasgn(this,s,b)
            s=s.subs;
            t=s{1};
            if ~ischar(t)
                t=mod(t,2);
                t=char(t+48);
            end
            set(this,t,b);
        end
        %% 加法
        function h=plus(f,g)
            if isnumeric(f)&&length(f)==1
                h=g;
                N=h.n;t=zeros(1,N);
                set(h,t,get(h,t)+f);
                return
            end
            if isnumeric(g)&&length(g)==1
                h=plus(g,f);
                return
            end
            if f.n~=g.n
                error('the addition are not on same group')
            end
            h=CZ_2nBig(f.n);
            Index=[cell2mat(keys(f.c).');cell2mat(keys(g.c).')];
            Index=unique(Index,'rows');
            for i=1:(size(Index,1))
                t=Index(i,:);
                if get(f,t)+get(g,t)~=0
                    set(h,t,get(f,t)+get(g,t));
                end
            end
        end
        
        function L=length(this)
            L=length(this.c);
        end
        
        %% 减法
        function h = minus(f,g)
            if isnumeric(f)&&length(f)==1
                h=minus(CZ_2nBig(g.n),g)+f;
                return
            end
            if isnumeric(g)&&length(g)==1
                h=plus(f,-g);
                return
            end
            if f.n~=g.n
                error('the addition are not on same group')
            end
            h=CZ_2nBig(f.n);
            Index=[cell2mat(keys(f.c).');cell2mat(keys(g.c).')];
            Index=unique(Index,'rows');
            for i=1:(size(Index,1))
                t=Index(i,:);
                if (get(f,t)-get(g,t))~=0
                    set(h,t,get(f,t)-get(g,t));
                end
            end
        end
        %% 乘法
        function h = mtimes(f,g)
            %C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
            if isnumeric(f)&&length(f)==1% is  scalar
                h=CZ_2nBig(g.n);
                Index=cell2mat(keys(g.c).');
                for i=1:(size(Index,1))
                    t=Index(i,:);
                    set(h,t,f*get(g,t));
                end
                return
            end
            if isnumeric(g)&&length(g)==1% is  scalar
                h=g*f;
                return
            end
            N=f.n;
            h=CZ_2nBig(N);
            Index1=cell2mat(keys(f.c).');
            Index2=cell2mat(keys(g.c).');
            for i=1:size(Index1,1)
                x=Index1(i,:);
                x=x-48;
                for j=1:size(Index2,1)
                    y=Index2(j,:);
                    y=y-48;
                    t=mod(x+y,2);
                    set(h,t,get(h,t)+get(f,x)*get(g,y));
                end
            end
        end
        %% 乘法2 .*
        function h=times(f,g)
            %C = TIMES(A,B) is called for the syntax 'A .* B'
            h=f*g;
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
            K=keys(f.c);
            for i=K
                ind=i{1};
                ind=ind-48;
                y=get(f,ind);
                term=y*prod(x(:).^(ind(:)));
                S=S+term;
            end
        end
        
        %% Find
        function [subs,vals] = find(f)
            subs=keys(f.c);
            subs=cell2mat(subs.');
            vals=values(f.c);
            vals=cell2mat(vals.');
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
        
        %% IFFT
        function V=CZifft(f)
            [A,B]=find(f);
            N=f.n;
            A=A-48;
            A=A;
            ProInd=(1:N)-1;
            ProInd=2.^(ProInd(:));
            Ind=A*ProInd+1;
            N=f.n;
            F=zeros(ones(1,N)*2);
            for i=1:size(A,1)
                F(Ind)=B;
            end
            V=ifftn(F)*2^N;
            V=reshape(flipud(V(:)),2*ones(1,f.n));
        end
    end
    
    
end