function [W,H]=scRCMF_c(data,K)
%gap_cluster Compute the optimatical decomposition of single cell matrix
%INPUT:
%data : the single cell matrix
%K    :  ptimatical cluster number by BIC priniciple
times=20;[m,~]=size(data);X=data';
WH=cell(times,2);LAM=BIC_lam(data,K);
R=rand(m,K);r=zeros(1,times);I=eye(m,m);
for i=1:times
    i,
    [W,H]=scRNMF_cluster(X,K,R,LAM);
    r(i) = mean(mean(abs(X-W*H)))+LAM*mean(mean(abs(I-R*H)));%
    WH{i,1}=W;WH{i,2}=H;
end
[~,ind]=min(r);
W=WH{ind,1};H=WH{ind,2};

function LAM=BIC_lam(data,K)
%gap_cluster Compute the optimatical cluster number in single cell matrix
[m,~]=size(data);
lam=[0.001,0.01,0.1,1,10,100,1000];
times=20;X=data';R=rand(m,K);r=zeros(1,length(lam));
I=eye(m,m);
for j=1:length(lam)
   % j
    sj=0;
    for i=1:times
        [W,H]=rnmf2(X,K,R,lam(j));
        errorx = mean(mean(abs(X-W*H)))+lam(j)*mean(mean(abs(I-R*H)));%
        sj= sj-2 * log(errorx) + 4* log(m);%4 is the number of parameter
    end
    r(j)=sj/times;
end
[~,ind]=min(r);
LAM=lam(ind);%