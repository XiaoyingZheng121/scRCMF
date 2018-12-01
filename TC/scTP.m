function U=scTP(X2,Co_module,tc,tk,Hp)
%X2: cell matrix
% -- Co_module: Cell clusters corresponding key genes.
% --  tc: Transition states among between clusters
%--  tk: cluster
%-- Hp: probability matrix
%  Output
[~,scores,~] = pca(X2');
transMatrix = scores(:,1:3);
X1=transMatrix;
tw=tc;%X2=X1';
tc_ex=X1(tw,:);
center=zeros(size(X1,2),2);
xc=corr(X2);
for i=1:length(tk)
    x1=Co_module{tk(i),2};ci=xc(x1,x1);
    [~,ic]=max(mean(ci));
    center(:,i)=X1(x1(ic),:);
     %center(:,i)=mean(X2(:,x1)');
end
Data=tc_ex;H1=Hp(tk,tw);
 for i=1:2
     H1(:,2)=H1(:,2)/sum(H1(:,2));
 end
U=Fuzzycm2(Data',center,H1);
