% compare kapp or kmax to kcat

load('kapp.mat');
load('kcat.mat');
load('kmax.mat');

% %remove reactions with promiscuous enzymes and isozymes
% load('GEM-yeast-split.mat');
% model = model_split;
% rxn_kept = cell(0,1);
% for i = 1:length(model.rxns)
%     if ~contains(model.grRules(i),'(') && ~ismember(model.grRules(i),'')
%         if strlength(model.rxns{i}) == 6
%             if sum(ismember(model.grRules,model.grRules(i))) == 1
%                 rxn_kept = [rxn_kept;model.rxns(i)];
%             end
%         else
%             if sum(ismember(model.grRules,model.grRules(i))) == 2
%                 rxn_kept = [rxn_kept;model.rxns(i)];
%             end
%         end
%     end
% end
% idxkapp = ismember(kapp.rxn,rxn_kept);
% kapp.protein = kapp.protein(idxkapp,:);
% kapp.rxn = kapp.rxn(idxkapp,:);
% kapp.values = kapp.values(idxkapp,:);
% idxkmax = ismember(kapp4.rxn,rxn_kept);
% kapp4.fc = kapp4.fc(idxkmax,:);
% kapp4.max = kapp4.max(idxkmax,:);
% kapp4.max_cond = kapp4.max_cond(idxkmax,:);
% kapp4.maxFlux = kapp4.maxFlux(idxkmax,:);
% kapp4.maxkmax = kapp4.maxkmax(idxkmax,:);
% kapp4.minFlux = kapp4.minFlux(idxkmax,:);
% kapp4.minkmax = kapp4.minkmax(idxkmax,:);
% kapp4.pFBA = kapp4.pFBA(idxkmax,:);
% kapp4.protein = kapp4.protein(idxkmax,:);
% kapp4.rxn = kapp4.rxn(idxkmax,:);


%% color

maincolor = [240,59,32]/255;
heatmaplow = [255,237,160]/255;

%% correlate kcat with kmax (rxns that have non-zero kapps at >= 4 conditions)

figure();
line([-2 6],[-2 6],'Color','k');
hold on;
box on;
rxns = intersect(kcat.rxn,kapp4.rxn);
[~,p] = ismember(rxns,kcat.rxn);
x_kcat = kcat.value(p);
[~,q] = ismember(rxns,kapp4.rxn);
y_kmax = kapp4.max(q);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-0.5,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1.2,['p' ' = ' num2str(round(PVAL,9))],'Color','black','FontSize',6,'FontName','Helvetica');
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('all data','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);

figure();
rxns = intersect(kcat.rxn,kapp4.rxn);
[~,p] = ismember(rxns,kcat.rxn);
x_kcat = kcat.value(p);
[~,q] = ismember(rxns,kapp4.rxn);
y_kmax = kapp4.max(q);
idx_heterexp = kcat.HeterExp(p) ~= 1;
x_kcat = x_kcat(idx_heterexp);
y_kmax = y_kmax(idx_heterexp);
rxns = rxns(idx_heterexp);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
line([-2 6],[-2 6],'Color','k');
hold on;
box on;
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-0.5,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1.2,['p' ' = ' num2str(round(PVAL,10))],'Color','black','FontSize',6,'FontName','Helvetica');
rxn_bigdev = rxns(abs(log10(x_kcat./y_kmax)) > 2);
for i = 1:length(rxn_bigdev)
    idxtmp = ismember(rxns,rxn_bigdev{i});
    text(log10(x_kcat(idxtmp)),log10(y_kmax(idxtmp)),strrep(rxn_bigdev{i},'_',''),'Color','black','FontSize',5,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
end
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('homologous expression','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[400 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);

figure();
rxns = intersect(kcat.rxn,kapp4.rxn);
[~,p] = ismember(rxns,kcat.rxn);
x_kcat = kcat.value(p);
[~,q] = ismember(rxns,kapp4.rxn);
y_kmax = kapp4.max(q);
idx_heterexp = kcat.HeterExp(p) == 1;
x_kcat = x_kcat(idx_heterexp);
y_kmax = y_kmax(idx_heterexp);
rxns = rxns(idx_heterexp);
load('PTMinfo.mat'); % check PTMs
PTMinfo.protein = PTMinfo.protein(ismember(PTMinfo.type,'Phosphorylation')|ismember(PTMinfo.type,'Ubiquitination'));
ptmtf = zeros(length(rxns),1);
ptmprotlist = cell(0,1);
for i = 1:length(rxns)
    prottmp = kapp4.protein(ismember(kapp4.rxn,rxns(i)));
    prottmp = strrep(prottmp,'( ','');
    prottmp = strrep(prottmp,' )','');
    prottmp = split(prottmp,' or ');
    ptmtf(i) = any(ismember(prottmp,PTMinfo.protein));
    ptmprotlist = [ptmprotlist;prottmp(ismember(prottmp,PTMinfo.protein))];
end
ptmprotlist = unique(ptmprotlist);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
line([-2 6],[-2 6],'Color','k');
hold on;
box on;
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-0.5,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1.2,['p' ' = ' num2str(round(PVAL,3))],'Color','black','FontSize',6,'FontName','Helvetica');
rxn_bigdev = rxns(abs(log10(x_kcat./y_kmax)) > 2);
for i = 1:length(rxn_bigdev)
    idxtmp = ismember(rxns,rxn_bigdev{i});
    text(log10(x_kcat(idxtmp)),log10(y_kmax(idxtmp)),strrep(rxn_bigdev{i},'_',''),'Color','black','FontSize',5,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
end
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('heterologous expression','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);


figure();
line([-2 6],[-2 6],'Color','k');
hold on;
rxns = intersect(kcat.rxn,kapp4.rxn);
[~,p] = ismember(rxns,kcat.rxn);
x_kcat = kcat.value(p);
[~,q] = ismember(rxns,kapp4.rxn);
y_kmax = kapp4.max(q);
y_low = kapp4.minkmax(q);
y_high = kapp4.maxkmax(q);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0);
text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-0.5,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1.2,['p' ' = ' num2str(round(PVAL,9))],'Color','black','FontSize',6,'FontName','Helvetica');
for i = 1:length(y_kmax)
    line([log10(x_kcat(i)),log10(x_kcat(i))],[log10(y_low(i)),log10(y_high(i))],'Color',maincolor,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_low(i)),log10(y_low(i))],'Color',maincolor,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_high(i)),log10(y_high(i))],'Color',maincolor,'LineWidth',0.5);
end
box on;
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('data with FVA as error bar','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[800 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);


figure();
line([-2 6],[-2 6],'Color','k');
hold on;
rxns = intersect(kcat.rxn,kapp4.rxn);
[~,p] = ismember(rxns,kcat.rxn);
x_kcat = kcat.value(p);
[~,q] = ismember(rxns,kapp4.rxn);
y_kmax = kapp4.max(q);
y_low = kapp4.minkmax(q);
y_high = kapp4.maxkmax(q);
idx = log10(y_high)-log10(y_low) < 0.5; % remove data with log10(high/low) ≥ 0.5
x_kcat = x_kcat(idx);
y_kmax = y_kmax(idx);
y_low = y_low(idx);
y_high = y_high(idx);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0);
text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-0.5,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1.2,['p' ' = ' num2str(round(PVAL,8))],'Color','black','FontSize',6,'FontName','Helvetica');
for i = 1:length(y_kmax)
    line([log10(x_kcat(i)),log10(x_kcat(i))],[log10(y_low(i)),log10(y_high(i))],'Color',maincolor,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_low(i)),log10(y_low(i))],'Color',maincolor,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_high(i)),log10(y_high(i))],'Color',maincolor,'LineWidth',0.5);
end
box on;
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('data without high flux variability','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[1000 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);

%% saturation with growth rate
[num,txt,~] = xlsread('ProteomicsFlux.xlsx','Flux');
exRxnList = txt(2:end,1);
exFluxes = num;
expList = txt(1,3:end);
clear num txt;
expGR = exFluxes(ismember(exRxnList,'r_2111'),:);

list_cond = kapp.condition;
list_gr = zeros(1,length(list_cond));
list_satur = zeros(length(kapp4.max),length(list_cond));
[~,idxmax] = ismember(kapp4.rxn,kapp.rxn);
for i = 1:length(list_cond)
    condtmp = kapp.condition(i);
    list_gr(i) = mean(expGR(contains(expList,condtmp)));
    list_satur(:,i) = kapp.values(idxmax,i)./kapp4.max;
end

figure();
for i = 1:length(list_gr)
    x = list_gr(i);
    saturtmp = list_satur(:,i);
    saturtmp = saturtmp(saturtmp~=0);
    y = mean(saturtmp);
    hold on;
    scatter(x,y,15,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
    
end
box on;
xlim([0 0.45]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Average kapp/kmax','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[100 100 100 100]);
set(gca,'position',[0.3 0.2 0.6 0.7]);

grcutoff = 0.20;
satur_low = list_satur(:,list_gr <= grcutoff);
satur_low = satur_low(:);
satur_low = satur_low(satur_low~=0);
satur_high = list_satur(:,list_gr > grcutoff);
satur_high = satur_high(:);
satur_high = satur_high(satur_high~=0);

lenlow = length(satur_low);
lenhigh = length(satur_high);
maxnum = max([lenlow,lenhigh]);
satur_low = [satur_low;nan(maxnum-length(satur_low),1)];
satur_high = [satur_high;nan(maxnum-length(satur_high),1)];
figure();
h = boxplot([satur_low satur_high],'Notch','on','Symbol','.','OutlierSize',6,'Widths',0.3,'Colors',maincolor);
set(h,{'linew'},{0.5});
xlim([0.5 2.5]);
xticks([1 2]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'XTickLabel',{['slow(N = ' num2str(lenlow) ')'],['fast(N = ' num2str(lenhigh) ')']});
set(gca,'FontSize',6,'FontName','Helvetica');
set(gca,'XColor','k');
set(gca,'YColor','k');
ylabel('Individual kapp/kmax','FontSize',7,'FontName','Helvetica','Color','k');
ranksump = ranksum(satur_low,satur_high);
title(['rank sum test p' ' < 1e' num2str(ceil(log10(ranksump)))],'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[200 100 100 100]);
set(gca,'position',[0.25 0.2 0.6 0.7]);

