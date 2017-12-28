%% Setup paths and load data
clear;
clc;
echo on;
addpath('Data');
addpath('Cluster');
addpath('TC');

% Select and read data matrix in format: Genes x Cells 
cell_data = xlsread('cell_HEE.xls','Sheet1');%load HEE_matrix.mat;
data = cell_data;
%data=X;
% Load all gene names as string
load gene_name.mat;
% m=size(data,1);n=size(data,2);
% data=cell_data;
% %% Optional step: preprocess data by selecting a subset of genes
data(all(var(data,0,2)<0.5,2),:)=[];
% 
%% Step 1: Run scRCMF to identify clusters and subpopulation composition
K= [];    % K is the number of clusters: can be specified by user, or
            % if not given (K = []), it will be inferred

 [Cell_module,kopt] = scRCMF_cluster(data,K);
% Output
%   -- Cell_module: Cell clusters corresponding key genes.
%   --  kopt: Number of clusters
%   --  pp: OI
%   --  w: left matrix
%   --  h: right matrix

%% Step 2: Run scRCMF to identify and analysis transition states 
 tc=scRCMF_TC(Cell_module,kopt);
 % Output
%  -- tc£ºtransition cell between  cell clusters

TE=scRCMF_TE(tc,cell_data,Cell_module,kopt);
 % Output
%  -- TE: entropy of transition states corresponding cluster label

pt=scRCMF_TP(tc,cell_data,Cell_module,kopt);
%  -- pt:  transition probability from  transition state to corresponding cell cluster

%% Step 3: select differnential gene 
gene_set=gene_pattern(cell_data,Cell_module,tc,kopt,gene_name);



