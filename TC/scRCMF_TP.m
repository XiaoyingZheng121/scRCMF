function U=scRCMF_TP(X,H,Co_module,tc_module)
 % Input
% -- X  and H    : matrix and 
% --  cell_data : single cell matrix
% -- Co_module: Cell clusters corresponding key genes.
% --  tc_module: Transition states among between clusters
%  Output
k= size(H,1);
m=size(H,2);
Hp=H;
for u=1:k 
    Hp(u,:)=Hp(u,:)/max(Hp(u,:));
end
for v=1:m
    Hp(:,v)=Hp(:,v)/sum(Hp(:,v));%
end

U=cell(1,size(tc_module,1));
for i=1:size(tc_module,1)
    tk=tc_module{i,1};tc=tc_module{i,2};
    U{i,1}=tk;U{i,2}=scTP(X,Co_module,tc,tk,Hp);
end