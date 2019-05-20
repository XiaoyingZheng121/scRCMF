%% Setup paths and load data
clear;
clc;
echo on;
addpath('Data');
addpath('Cluster');
addpath('TC');

%input data
cell_data=load('sim_data.mat');

%group2=idx;
data = cell_data.X1;
%data=X;
% Load all gene names as string
%load gene_name.mat;

 %% Optional step: preprocess data by selecting a subset of genes
%data(all(var(data,0,2)<0.5,2),:)=[];
% 
%% Step 1: Run scRCMF to identify  low rank structure 
K=gap_cluster(data');    % K is the number of clusters: can be specified by user, or
            % if not given (K = []), it will be inferred

[W,H] = scRCMF_c(data,K);

% Output
%   --  kopt: Number of clusters
%   --  W: left matrix
%   --  H: right matrix

%% Step 2: Run scRCMF to identify subpopulations and transition states 
c=0.6;       % c is the thresold: can be specified by user
Co_module  = Cell_cluster(H,W,c);%
tc_module = identify_tc(H,c);
 % Output
%   -- Co_module: Cell clusters corresponding key genes.
%  -- tc_module£:transition states between  cell clusters

%% Step 3: Run scRCMF to analyse transition states 
TE=scRCMF_TE(H);
 % Output
%  -- scEntropy: Entropy of transition states corresponding cluster label

U=scRCMF_TP(data,H,Co_module,tc_module);
%  -- scTP: transition probability from  transition state to corresponding cell cluster

%% Step 4: select differnential gene 
%gene_set=gene_pattern(cell_data,Co_module,tc_module,K,gene_name);%HEE data



