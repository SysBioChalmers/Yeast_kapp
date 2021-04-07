% second max kapp vs kcat
load('kapp.mat');

kapp4.max = zeros(0,1);
kapp4.rxn = cell(0,1);
kapp4.protein = cell(0,1);
kapp4.fc = zeros(0,1);
for i = 1:length(kapp.rxn)
    values_tmp = kapp.values(i,:);
    if sum(values_tmp > 0) >= 4
        kapp4.rxn = [kapp4.rxn;kapp.rxn(i)];
        kapp4.protein = [kapp4.protein;kapp.protein(i)];
        values_tmp(values_tmp == max(values_tmp)) = 0;
        kapp4.max = [kapp4.max;max(values_tmp)];
        kapp4.fc = [kapp4.fc;max(values_tmp)/min(values_tmp)];
    end
end

load('kcat.mat');

%% color
% maincolor = [124,81,161]/255;
% heatmaplow = [218,218,235]/255;

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
text(2,-1.2,['p' ' = ' num2str(round(PVAL,10))],'Color','black','FontSize',6,'FontName','Helvetica');
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('all data','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 120 120]);
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
rxn_bigdev = rxns(abs(log10(x_kcat./y_kmax)) > quantile((abs(log10(x_kcat./y_kmax))),0.95));
for i = 1:length(rxn_bigdev)
    idxtmp = ismember(rxns,rxn_bigdev{i});
    text(log10(x_kcat(idxtmp)),log10(y_kmax(idxtmp)),strrep(rxn_bigdev{i},'_',''),'Color','black','FontSize',5,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
end
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('data without heterologous expression','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[400 200 120 120]);
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
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
line([-2 6],[-2 6],'Color','k');
hold on;
box on;
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-0.5,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1.2,['p' ' = ' num2str(round(PVAL,10))],'Color','black','FontSize',6,'FontName','Helvetica');
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('data of heterologous expression','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 200 120 120]);
set(gca,'position',[0.2 0.2 0.7 0.7]);
