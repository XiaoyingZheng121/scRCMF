function K=gap_cluster(data)
%gap_cluster Compute the optimatical cluster number in single cell matrix
%data is the single cell matrix
% PCA and visualization
[~,score] = pca(zscore(data)');
score=score(:,1:2);
eva = evalclusters(score,@kmedoids,'gap','KList',[1:10],'Distance','correlation');
K = eva.OptimalK;
