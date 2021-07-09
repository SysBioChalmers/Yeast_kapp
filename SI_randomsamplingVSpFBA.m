% correlate random sampling-based kmax with pFBA-based kmax

load('kmax.mat');
rxn1 = kapp4.rxn;
kmax1 = kapp4.max;
clear kapp4;


CalType = 'median';
% CalType = 'mean';

load(['kmax_sampling_' CalType '.mat']);
rxn2 = kapp4.rxn;
kmax2 = kapp4.max;
clear kapp4;

if strcmp(CalType,'median')
    maincolor = [215,48,39]/255;
elseif strcmp(CalType,'mean')
    maincolor = [69,117,180]/255;
end

figure();
line([-3 6],[-3 6],'Color','k');
hold on;
box on;
rxns = intersect(rxn1,rxn2);
[~,p] = ismember(rxns,rxn1);
x_kcat = kmax1(p);
[~,q] = ismember(rxns,rxn2);
y_kmax = kmax2(q);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
scatter(log10(x_kcat),log10(y_kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
text(2,0,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-1,['N = ' num2str(length(x_kcat))],'Color','black','FontSize',6,'FontName','Helvetica');
text(2,-2,['p' ' = ' num2str(round(PVAL,9))],'Color','black','FontSize',6,'FontName','Helvetica');
xlim([-3 6]);
ylim([-3 6]);
xticks(-3:3:6);
yticks(-3:3:6);
% title('all data','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kmax (pFBA) (/s)','FontSize',7,'FontName','Helvetica');
ylabel(['log10 kmax (' CalType ' sampling) (/s)'],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 110 110]);
set(gca,'position',[0.2 0.2 0.7 0.7]);


