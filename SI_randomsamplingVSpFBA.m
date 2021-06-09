% correlate random sampling-based kmax with pFBA-based kmax

load('kmax.mat');
rxn1 = kapp4.rxn;
kmax1 = kapp4.max;
clear kapp4;

load('kmax_sampling_median.mat');
rxn2 = kapp4.rxn;
kmax2 = kapp4.max;
clear kapp4;

maincolor = [240,59,32]/255;

figure();
line([-2 6],[-2 6],'Color','k');
hold on;
box on;
rxns = intersect(rxn1,rxn2);
[~,p] = ismember(rxns,rxn1);
x_kcat = kmax1(p);
[~,q] = ismember(rxns,rxn2);
y_kmax = kmax2(q);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-0.5,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1.2,['p' ' = ' num2str(round(PVAL,9))],'Color','black','FontSize',6,'FontName','Helvetica');
% rxn_bigdev = rxns(abs(log10(x_kcat./y_kmax)) > quantile((abs(log10(x_kcat./y_kmax))),0.9));
% for i = 1:length(rxn_bigdev)
%     idxtmp = ismember(rxns,rxn_bigdev{i});
%     text(log10(x_kcat(idxtmp)),log10(y_kmax(idxtmp)),strrep(rxn_bigdev{i},'_',''),'Color','black','FontSize',5,'FontName','Helvetica','HorizontalAlignment','center','VerticalAlignment','bottom');
% end
xlim([-2 5]);
ylim([-2 5]);
xticks(-2:1:5);
yticks(-2:1:5);
title('all data','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kmax (pFBA) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (random sampling) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);


