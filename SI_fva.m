% compare kapp or kmax to kcat

load('kapp.mat');
load('kcat.mat');
load('kmax.mat');

%% color

maincolor = [240,59,32]/255;
heatmaplow = [255,237,160]/255;

%% correlate kcat with kmax (rxns that have non-zero kapps at >= 4 conditions)

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
idx = log10(y_high)-log10(y_low) < 1; % remove data with log10(high/low) ≥ 1
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
idx_heterexp = kcat.HeterExp(p) == 1;
x_kcat = x_kcat(idx_heterexp);
y_kmax = y_kmax(idx_heterexp);
rxns = rxns(idx_heterexp);
y_low = y_low(idx_heterexp);
y_high = y_high(idx_heterexp);
idx = log10(y_high)-log10(y_low) < 1; % remove data with log10(high/low) ≥ 1
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
title('heterologous data without high flux variability','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[1200 200 110 110]);
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
idx_heterexp = kcat.HeterExp(p) ~= 1;
x_kcat = x_kcat(idx_heterexp);
y_kmax = y_kmax(idx_heterexp);
rxns = rxns(idx_heterexp);
y_low = y_low(idx_heterexp);
y_high = y_high(idx_heterexp);
idx = log10(y_high)-log10(y_low) < 1; % remove data with log10(high/low) ≥ 1
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
title('homologous data without high flux variability','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[1400 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);