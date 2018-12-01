function [ Co_module ] = identify_tc(Hp,c)
% Compute the transition states among cell clusters 
 % Input
%  -- Hp: probability matrix of cell
%  -- c£ºthreshold

%
% Output rule 2


K= size(Hp,1);
m=size(Hp,2);
for u=1:K 
    Hp(u,:)=Hp(u,:)/max(Hp(u,:));
end
for v=1:m
    Hp(:,v)=Hp(:,v)/sum(Hp(:,v));%
end

H1=Hp;
H1(H1>c)=1;

Co_module={};
t=0;hi=zeros(K,m);
for v=1:m  
    lv=H1(:,v);lf=find(lv==1);
    if isempty(lf)
    [~,lv]=sort(lv,'descend');
   hi(lv(1:2),v)=1; 
    end
end

for i=1:(K-1)
    for j=(i+1):K
        h=hi(i,:)+hi(j,:);
        % h=abs(Hij);%w=var(Wij);
        %m1=mean(h);v1=std(h);%m2=mean(h);v2=std(h);
%         r1=find(hi(i,:)==1); 
%         r2=find(hi(j,:)==1);       
        r1=find(h==2);        
        %Co_module{t,1}=c1;
        if ~isempty(r1)
                 t=t+1;
            Co_module{t,2}=r1;Co_module{t,1}=[i,j];
            %Co_module{i,3}=c2'; Co_module{i,4}=c3';       
        end
        %    Subpattern1{i}=X1(r,c1');
        %    Subpattern2{i}=X2(r,c2');
        %    Subpattern3{i}=X2(r,c3');
    end
end

