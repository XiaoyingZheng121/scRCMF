function [W,H]=rnmf2(X,K,R,l)
%Solving the new optimization function added to another random constrain about right matrix

%    min_{W,H}  ||X-W*H|| + ||I-R*H||
%    s.t.       W>=0; 
%               H>=0; 
%              
% INPUT:
% X  (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K        : Number of components
% R        : constraint of decompose right matrix

%
% OUTPUT:
% W        : N x K matrix
% H        : K x M matrix
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxiter=1000;    %Maximum number of iterations to run
speak=0;         %prints iteration count and changes in connectivity matrix
print_iter = 50; % iterations between print on screen and convergence test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%l=1;
if min(min(X)) < 0
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W and H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n,m]=size(X);
W=rand(n,K);
H=rand(K,m);
T=rand(m,m);
%r2=rand(K,m);
% use W*H to test for convergence
Xr_old=W*H;I=eye(m);
for iter=1:maxiter
    % Euclidean multiplicative method
       %T=T.*(2*R*H)./(I+R*(H*H')*R');
    H = H.*(W'*X+l*R')./((W'*W*H+l*(R'*R)*H)+eps);%a
    %H = H.*(W'*X+R)./((W'*W*H+R*R'*H)+eps);%add to the constraint
    W = W.*(X*H')./(W*(H*H')+eps);
    % H = H.*(W'*X+R)./((W'*W*H+R*R'*H)+eps);
    % print to screen
    if (rem(iter,print_iter)==0) && speak,
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist = nmf_euclidean_dist(X,W*H);
        % errorx = mean(mean(abs(X-W*H)))/mean(mean(X));
        errorx = mean(mean(abs(X-W*H)))/mean(mean(X))+mean(mean(abs(I-R*H)))/mean(mean(X));%
        disp(['Iter = ',int2str(iter),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end

function err = nmf_euclidean_dist(X,Y)

err = sum(sum((X-Y).^2));
