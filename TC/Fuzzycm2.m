function U=Fuzzycm2(Data,P0,Hp)

% 
% Input: 
% --  Data  : expression of transition cell, row is the label cell
% --  P0    : cluster center
% Output: 
% --  U     : transition probability from  transition state to corresponding cell cluster
v=P0;tc_ex=Data;
N=size(tc_ex,2);%sample number
k=size(v,2);%
p0=0;iter=0;
for i=1:k
    for j=1:N
        pi=Hp(i,j).^2*norm(v(:,i)-tc_ex(:,j),2)^2;
        p0=p0+pi;
    end
end

while true  
    iter=iter+1; 
     for i=1:k
         ti=zeros(length(v(:,i)),1);
         for j=1:N
             ti=ti+Hp(i,j).^2*tc_ex(:,j);
         end
     v(:,i)=ti./(sum(Hp(i,:).^2)+eps);
     end
    % compute U
    flag=0;
    for i=1:k
        for j=1:N
            if norm(v(:,i)-tc_ex(:,j))<0.01            
                Hp(:,j)=zeros(1,k);
                Hp(i,j)=1;flag=1;
            end
        end
    end
  if flag==0
      for i=1:k
          for j=1:N
              xk=tc_ex(:,j);
              t1=norm(v(:,i)-xk,2)^2;t2=0;
             % xx=repmat(xk,1,3);
             for l=1:k
                 t2=t2+t1/norm(v(:,l)-xk,2)^2;
             end
             % t2=norm(v-xx,'fro')^2;
          Hp(i,j)=1.0./t2;
          end
      end
  end
  P=0;
for i=1:k
    for j=1:N
        pi=Hp(i,j).^2*norm(v(:,i)-tc_ex(:,j),2)^2;
        P=P+pi;
    end
end
    % stop criterion
     if  abs(P-p0)<0.1       break
     end
%        if  abs(Obj_Fcn(iter)-Obj_Fcn(iter-1)) <epsm      break
%      end
    p0=P;
end
U=Hp;