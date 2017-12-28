function overlap_cell=scRCMF_TC(Cell_module,k)
% Input
%   -- Cell_module: Cell clusters and gene clusters.
%   --  k: Number of clusters

%
% Output
% -- overlap cells£ºtransition cell between  cell clusters
overlap_cell=[];

for i=1:(k-1)
    v1=Cell_module{i,2}';
    for j=(i+1):k
        if j~=i
            v2=Cell_module{j,2}';
            if ~isempty(intersect(v1,v2)) && ~isempty(v1) && ~isempty(v2)
                overlap_cell=[overlap_cell,intersect(v1,v2)];%setdiff
                %intersect(v1,v2),
            end
        end
    end
end
overlap_cell=unique(overlap_cell);


