function entropy = scRCMF_TE(H)
% % compare the entropy of single cell

k= size(H,1);
m=size(H,2);
for u=1:k 
    H(u,:)=H(u,:)/max(H(u,:));
end

for v=1:m
    H(:,v)=H(:,v)/sum(H(:,v));%
end

entropy=zeros(1,m);
for i=1:m
    entropy(i)=-sum(H(:,i).*log(H(:,i)+eps));
end