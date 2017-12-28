function [kopt,p]=identify_k(X,nk)
% This is the  function to identify the optimal number of cluster with AEC index
%
% Input:
%        X: A m*n matrix with m rows (genes) and n columns (cells).
%       nk: max value of Number of cluster.
%
% Output:
%           kopt: optimal number of cluster with minimum of AEC index

%
    if isempty(nk)
        nk =20;
    end
    
p=zeros(nk,1);

tt0=2;tt1=2;%threshold of gene and cell
n=size(X,2);
for k=2:nk
r=rand(k,n);
    [W, H]=rnmf(X,k,r); %nonnegative matrix factorization
    [ Co_module ] = Classification(W, H, tt0, tt1);
    cell_s=[];
    for i=1:k
        v1=Co_module{i,2}';% cell subpopulation
         st=[];x1=Co_module{i,2};st=[st;x1];
        for j=(i+1):k
            if j~=i
                v2=Co_module{j,2}';
                if ~isempty(setdiff(v1,v2))
                    cell_s=[cell_s,intersect(v1,v2)];
                end
            end
        end
    end
    %if length(unique(cell_s))>3
    p(k)=length(cell_s)/nchoosek(k,2);%average coefficient of variation
   % end
   % if p(k)>0
    fprintf('%d, %8.6f\n',k,p(k));
   % end
end
%find the optimal number of cluster
p(p(:)==0)=inf;
[~,kopt]=min(p(:));




