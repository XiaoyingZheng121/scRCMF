function geneset=gene_pattern(cell_data,Cell_module,tc,kopt,gene_name)
 % Input
% -- tc         : transition cell between  cell clusters
% -- cell_data  : single cell matrix
% -- Cell_module: Cell clusters corresponding key genes.
% -- kopt       : Number of clusters
%-- gene_name   : All gene names
m=size(cell_data,1);n=size(cell_data,2);
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
lm=unique(sort(lt,2),'rows');lo=setdiff(1:size(Cell_module,1),lm);
% key gene in first cluster pattern
gt1=0;tq1=0;
% x2=Cell_module{lm(1),2};x2=setdiff(x2,tc);
x3=Cell_module{lm(2),2};x3=setdiff(x3,tc);

%adjust parameter
r1=0.5;%number percent of key gene expression
r2=0.5;%variance of key gene expression
r3=2;  %control mean gene expression of specific cell cluster
r4=1;%control mean gene expression of other cell cluster
for i=1:m
    j=1;
    ee=i;   le=find(cell_data(ee,:)>0.1);
    x1=Cell_module{j,2};x1=setdiff(x1,tc);
    l1=mean(cell_data(ee,x1));e1=var(cell_data(ee,x1));
    l3=mean(cell_data(ee,x3)); %l2=mean(cell_data(ee,x2));
    if length(le)/n>r1&&l1>r3&&e1<r2&&l3<r4
        gt1=[gt1,i];tq1=[tq1,mean(cell_data(ee,tc))];
    end
end

gt1(1)=[];tq1(1)=[];
[~,ic] = sort(tq1,'descend');
gene1=gt1(ic);
% key gene in second cluster pattern

gt1=0;tq1=0;
x2=Cell_module{lo,2};x2=setdiff(x2,tc);
%x3=Cell_module{lm(2),2};x3=setdiff(x3,tc);
for i=1:m
    j=2;
    ee=i;   le=find(cell_data(ee,:)>0.1);
    x1=Cell_module{j,2};x1=setdiff(x1,tc);
    l1=mean(cell_data(ee,x1));e1=var(cell_data(ee,x1));
    l2=mean(cell_data(ee,x2));%l3=mean(cell_data(ee,x3));
    if length(le)/n>r1&&l1>r3&&e1<r2&&l2<r4
        gt1=[gt1,i];tq1=[tq1,mean(cell_data(ee,tc))];
    end
end
gt1(1)=[];tq1(1)=[];
[~,ic] = sort(tq1,'descend');
gene2=gt1(ic);

% key gene in third cluster pattern
gt1=0;tq1=0;
x2=Cell_module{lm(1),2};x2=setdiff(x2,tc);
x3=Cell_module{lo,2};x3=setdiff(x3,tc);
for i=1:m
    j=3;
    ee=i;   le=find(cell_data(ee,:)>0.1);
    x1=Cell_module{j,2};x1=setdiff(x1,tc);
    l1=mean(cell_data(ee,x1));e1=var(cell_data(ee,x1));
   l3=mean(cell_data(ee,x3));% l2=mean(cell_data(ee,x2));
    if length(le)/n>r1&&l1>r3&&e1<r2&&l3<r4% &&l3<1 %mean(cell_data(ee,tc))>3var(cell_data(ee,tc))<0.1%ww(i,3)==max(ww(i,:))  &&
        gt1=[gt1,i];tq1=[tq1,mean(cell_data(ee,tc))];
    end
end
gt1(1)=[];tq1(1)=[];
[~,ic] = sort(tq1,'descend');

gene3=gt1(ic);


%show key gene for each pattern
gene=[gene1,gene2,gene3];gene_name=upper(gene_name);
gene=unique(gene);
%geneset=gene_name(gene);
  l=30;
 %gene=[gene1(1:ll),gene2(1:ll),gene3(1:ll)];
 geneset=gene_name(gene(1:l));
% 
% dataForHeatmap = y;idxMarkerTop10=upper(geneset);
% dataForHeatmapZscore = zscore(dataForHeatmap,[],2);
% imagesc(dataForHeatmapZscore);
% %imagesc(y);
% colormap(redbluecmap)
% set(gca,'Ytick',1:length(idxMarkerTop10))
% set(gca,'YtickLabel',geneset,'FontName','Arial','FontSize',8)
