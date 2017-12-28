function [U,P,Dist,Cluster_Res,Obj_Fcn,iter]=Fuzzycm(Data,P0)
% 模糊 C 均值聚类 FCM: 从指定初始聚类中心开始迭代
% 
% Input: 
% --  Data  : expression of transition cell, row is the label cell
% --  P0    : cluster center
% Output: 
% --  U     : transition probability from  transition state to corresponding cell cluster

if  nargin<5
    epsm=1.0e-6;  
end
if  nargin<4
    M=2;
end
if  nargin<3
    plotflag=0;
end

[N,S] = size(Data);
m = 2/(M-1);
iter = 0;
C=size(P0,1);
Dist(C,N)=0;
U(C,N)=0;
P(C,S)=0;

while true  
    iter=iter+1;  
    % compute U
    for i=1:C
        for j=1:N
            Dist(i,j)=fuzzydist(P0(i,:),Data(j,:));
        end
    end         
    U=1./(Dist.^m.*(ones(C,1)*sum(Dist.^(-m))));      
    % update P
    Um=U.^M;
    P=Um*Data./(ones(S,1)*sum(Um'))';   
    % 
    if  nargout>4 | plotflag
        Obj_Fcn(iter)=sum(sum(Um.*Dist.^2));
    end
    % stop criterion
     if  norm(P-P0,Inf)        break
     end
%        if  abs(Obj_Fcn(iter)-Obj_Fcn(iter-1)) <epsm      break
%      end
    P0=P;
end

% result
if  nargout > 3
    res = maxrowf(U);
    for c = 1:C
        v = find(res==c);
        Cluster_Res(c,1:length(v))=v;
    end
end
% plot
% if  plotflag
%     fcmplot(Data,U,P,Obj_Fcn);
% end