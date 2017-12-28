function pt=scRCMF_TP(tc,cell_data,Cell_module,kopt)
 % Input
% -- tc         : transition cell between  cell clusters
% --  cell_data : single cell matrix
% -- Cell_module: Cell clusters corresponding key genes.
% --        kopt: Number of clusters
%  Output

% --          pt:  transition probability from  transition state to cell cluster
m=size(cell_data,1);n=size(cell_data,2);
for i=1:n
    cell_data(:,i)=cell_data(:,i)./sum(cell_data(:,i));
end

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
c2=Cell_module{lm(1),2};c3=Cell_module{lm(2),2};
C2=setdiff(c2,tc);C3=setdiff(c3,tc);
pt=zeros(length(tc),2); xc=[C2;C3];
xc(1:length(C2))=1;xc(length(C2)+1:end)=2;
cx=unique(xc);
% xc=cell_state([x1;y1]);xx=unique(xc);%if know the stage

%fuzzy membership analysis
center=zeros(length(cx),m);
for j=1:length(cx)
    x1= xc==cx(j);x_c=mean(cell_data(:,x1),2);
    center(j,:)=x_c;
end

center=zeros(length(cx),m);
for j=1:length(cx)
    x1= xc==cx(j);x_c=mean(cell_data(:,x1),2);
    center(j,:)=x_c;
end
tc_ex=cell_data(:,tc)';
U=Fuzzycm(tc_ex,center);
U=U';
pt=U;
% plot
% lgp={'C2','C3'};
% bar(U)
% set(gca,'Xtick',1:length(tc));
% set(gca,'XTickLabel',tc)%HEE
% legend(lgp,'FontSize',10,'Location','best');%,'Orientation','horizontal');
% set(gca,'FontName','Times New Roman');
% set(gca,'FontSize',10);
% ylabel('p');
