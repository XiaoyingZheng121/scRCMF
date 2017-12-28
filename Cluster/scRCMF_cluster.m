function [Cell_module,kopt] = scRCMF_cluster(X,k)
% RNMF identifies clusters and subpopulations from single cell data
%
% Input
%   --   X:
%       a m*n single cell data matrix with m rows(genes) and n columns(cells)
%   --   k:
%       Number of Cluster specified by user, if k = [] (default), then the algorithm will compute the number.

%
% Output
%   --  Cell_module: cell matrix with k subpopulations corresponding key genes
%   --  kopt     : optimal k with the AEC index

if nargin==2 && ~isempty(k)
    kopt=k;pp=100;
end
if isempty(k)  ||  nargin==1
    nk=[];
    [kopt,pp]=identify_k(X,nk);
end
eps=0.1;q=pp;
tt0=2;tt1=2;%threshold of gene and cell
n=size(X,2);st=0;
c1=0.8;c2=1.2;pp=min(pp)+eps;p=pp+eps;
kopt,
while p>pp  ||isempty(cell_s)|| st<c1*n || st>=c2*n
    R=rand(kopt,n);
    st,
    [W, H]=rnmf(X,kopt,R);
    Cell_module = Classification(W, H, tt0, tt1);
    st=[];
    for i=1:kopt
        x1=Cell_module{i,2};st=[st;x1];
    end
    st=length(unique(st));
        cell_s=[];
    for i=1:kopt
        v1=Cell_module{i,2}';% cell subpopulation
        for j=(i+1):kopt
            if j~=i
                v2=Cell_module{j,2}';
                if ~isempty(setdiff(v1,v2))
                    cell_s=[cell_s,intersect(v1,v2)];
                end
            end
        end
    end
    p=length(cell_s)/nchoosek(kopt,2);
end
q(q(:)==inf)=0;





