function TE=scRCMF_TE(tc,cell_data,Cell_module,kopt)
% Input
% -- tc         : transition cell between  cell clusters
% --  cell_data : single cell matrix
% -- Cell_module: Cell clusters corresponding key genes.
% --        kopt: Number of clusters
%  Output

% --       TE: entropy of transition states corresponding cluster label
lt=zeros(length(tc),2);
for i=1:length(tc)
    for j=1:kopt
        x2=Cell_module{j,2};
        if ismember(tc(i),x2)
            l=j;break
        end
    end
    lt(i,1)=l;
    for u=(l+1):kopt
        x2=Cell_module{u,2};
        if ismember(tc(i),x2)
            l=u;break
        end
    end
    lt(i,2)=l;
end
lm=unique(sort(lt,2),'rows');

ii=0;icc=0;
ex= mapminmax(cell_data',0,1);ex=ex';%normalization
cell_ex=ex;eps=0.001;
for i=1:kopt
    x2=Cell_module{i,2}; x1=Cell_module{i,1}; id=0;
    h=1;
    for j=1:length(x2)
        if ~ismember(x2(j),tc)
            ce=cell_ex(x1,x2(j));
            %d1=var(ce);  l1=hist(ce);
            % d2=-sum((l1/sum(l1)+eps).*log2(l1/sum(l1)+eps));
            d3=-sum((ce+eps).*log2(ce+eps));   h=h+1;
            id=[id;d3];
        end
    end
    id(1,:)=[];
    ii=[ii;id];icc=[icc;i+zeros(h-1,1)];
end
ii(1,:)=[];icc(1)=[];
ss=zeros(length(tc),1);
for i=1:length(tc)
    t=1;s=0;
    for j=1:kopt
        x1=Cell_module{j,2};
        if ismember(tc(i),x1)
            l=j;break;
        end
    end
    xe=Cell_module{l,1};  %xe=dl(xe);
    x1=Cell_module{l,2};
    ce=cell_ex(xe,tc(i));  l1=hist(ce);
    %s(1)=var(ce);s(2)=-sum((l1/sum(l1)+eps).*log2(l1/sum(l1)+eps));
    s=-sum((ce+eps).*log2(ce+eps));
    for u=(l+1):kopt
        x2=Cell_module{u,1};   %  x2=dl(x2);
        xy=Cell_module{u,2};
        if ismember(tc(i),x1) && ismember(tc(i),xy)%找到不同的类
            t=t+1; ce=cell_ex(xe,tc(i)); l2=hist(ce);po=u;
            % s(1)=s(1)+var(ce);s(2)=s(2)-sum((l2/sum(l2)+eps).*log2(l2/sum(l2)+eps));
            s=s-sum((ce+eps).*log2(ce+eps));
        end
    end
    ss(i)=s/t;% average entropy
    for tt=1:size(lm,1)
        if l==min(lm(tt,:))&& po==max(lm(tt,:))
            icc=[icc;mean(lm(tt,:))+kopt];
        end
    end
    
end
ii=[ii;ss];
TE=[ii,icc];
%plot entropy of cells 
%boxplot(ii,icc);
colorTSNE={'r','b','g','k'};
dataMarkersTSNE=ii;groupFocus=icc;
boxplot(dataMarkersTSNE,groupFocus);
h1 = findobj(gca,'tag','Median');
set(h1,{'linew'},{2.5})
set(h1,'Color','k')
h2 = findobj(gca,'Tag','Box');
for j=1:length(h2)
    patch(get(h2(j),'XData'),get(h2(j),'YData'),colorTSNE{j},'FaceAlpha',0.5);
end
ylabel('scEntropy');
hold on
box on;
set(gca,'LineWidth',1.0);
set(gca,'XTickLabel',{'C1','C2','C3','TC'})%HEE
