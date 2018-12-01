function cell_module=Cell_cluster(H,W,c)
% max_cluster identifies modelu from single cell data
%
% Input
%   --   H and W:
%       low lank structure of single cell data matrix 
%
% Output
%   -- cell_module: cell matrix with cell subpopulations corresponding key genes


n=size(W,1);k= size(H,1);
m=size(H,2);z1=zeros(k,m);
for u=1:k 
    H(u,:)=H(u,:)/max(H(u,:));
end

for v=1:m
    H(:,v)=H(:,v)/sum(H(:,v));%
end

for i=1:m
    Hi=find(H(:,i)>c);
    if ~isempty(Hi)
        [~,li]=sort(H(Hi,i));
        z1(Hi(li(end)),i)=1;
    end
end
z2=zeros(n,k);
% for i=1:m
%     [~,j]=max(H(:,i));
%     z1(j,i)=1;
% end
for i=1:n
    [~,j]=max(W(i,:));
    z2(i,j)=1;
end
for i=1:k
    cell_module{i,2}=find(z1(i,:)==1);
    cell_module{i,1}=find(z2(:,i)==1);
end

