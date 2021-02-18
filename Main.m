load('kapp.mat');
load('kcat.mat');
load('kmax.mat');


%% correlate kcat with condition-specific kapp
figure();
for i = 1:length(kapp.condition)
    kapptmp = kapp.values(:,i);
    idxtmp = kapptmp~=0;
    kapprxn = kapp.rxn(idxtmp);
    kapplist = kapptmp(idxtmp);
    rxnstmp = intersect(kcat.rxn,kapprxn);
    [~,p] = ismember(rxnstmp,kcat.rxn);
    x_kcat = kcat.value(p);
    [~,q] = ismember(rxnstmp,kapprxn);
    y_kapp = kapplist(q);
    [RHOtmp,PVALtmp] = corr(log10(x_kcat),log10(y_kapp),'Type','Pearson');
    subplot(5,7,i);
    line([-4 6],[-4 6],'Color','k');
    hold on;
    box on;
    scatter(log10(x_kcat),log10(y_kapp),3,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
    xlim([-4 6]);
    ylim([-4 6]);
    xticks([]);
    yticks([]);
	text(-3.5,4.7,num2str(i),'FontSize',6,'FontName','Helvetica','FontWeight','bold');
    text(-2,-2.5,['R^2' '= ' num2str(round(RHOtmp^2,2))],'FontSize',6,'FontName','Helvetica');
end
set(gcf,'position',[400 400 360 260]);

%% correlate kapps between conditions
rho2data = zeros(length(kapp.condition),length(kapp.condition));
pdata = zeros(length(kapp.condition),length(kapp.condition));
for i = 1:length(kapp.condition)
    kapp_i = kapp.values(:,i);
    idx_i = kapp_i~=0;
    rxn_i = kapp.rxn(idx_i);
    kapp_i = kapp_i(idx_i);
    for j = 1:length(kapp.condition)
        kapp_j = kapp.values(:,j);
        idx_j = kapp_j~=0;
        rxn_j = kapp.rxn(idx_j);
        kapp_j = kapp_j(idx_j);
        rxn_ij = intersect(rxn_i,rxn_j);
        [~,p] = ismember(rxn_ij,rxn_i);
        idata = kapp_i(p);
        [~,q] = ismember(rxn_ij,rxn_j);
        jdata = kapp_j(q);
        [RHOtmp,PVALtmp] = corr(log10(idata),log10(jdata),'Type','Pearson');
        rho2data(i,j) = RHOtmp^2;
        pdata(i,j) = PVALtmp;
    end
end

max(max(pdata));

figure();
xlbl = num2cell(1:1:length(kapp.condition));
xlbl_tmp = cellfun(@(x) [num2str(x),'.'],xlbl,'UniformOutput',false);
ylbl = strcat(xlbl_tmp,kapp.condition);
minclr = [189,215,231]/255;
maxclr = [8,81,156]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(xlbl,ylbl,rho2data,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','k');
title(['R^2 between kapps of ' num2str(length(kapp.condition)) ' conditions']);
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[820 320 330 330]);
set(gca,'position',[0.25 0.2 0.55 0.55]);

%% correlate kcat with kmax (rxns that have non-zero kapps at >= 3 conditions)

figure();
line([-2 6],[-2 6],'Color','k');
hold on;
rxns = intersect(kcat.rxn,kapp4.rxn);
[~,p] = ismember(rxns,kcat.rxn);
x_kcat = kcat.value(p);
[~,q] = ismember(rxns,kapp4.rxn);
y_kmax = kapp4.max(q);
[RHO,PVAL] = corr(log10(x_kcat),log10(y_kmax),'Type','Pearson');
scatter(log10(x_kcat),log10(y_kmax),15,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
box on;
xlim([-2 6]);
ylim([-2 6]);
xticks([-2 0 2 4 6]);
yticks([-2 0 2 4 6]);
title(['R^2' '= ' num2str(round(RHO^2,3)) '; p' ' = ' num2str(round(PVAL,2)) '; N' ' = ' num2str(length(x_kcat))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 120 120]);
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
scatter(log10(x_kcat),log10(y_kmax),15,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0);
for i = 1:length(y_kmax)
    line([log10(x_kcat(i)),log10(x_kcat(i))],[log10(y_low(i)),log10(y_high(i))],'Color',[8,81,156]/255,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_low(i)),log10(y_low(i))],'Color',[8,81,156]/255,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_high(i)),log10(y_high(i))],'Color',[8,81,156]/255,'LineWidth',0.5);
end
box on;
xlim([-2 6]);
ylim([-2 6]);
xticks([-2 0 2 4 6]);
yticks([-2 0 2 4 6]);
title(['R^2' '= ' num2str(round(RHO^2,3)) '; p' ' = ' num2str(round(PVAL,2)) '; N' ' = ' num2str(length(x_kcat))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[500 200 120 120]);
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
scatter(log10(x_kcat),log10(y_kmax),15,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0);
for i = 1:length(y_kmax)
    line([log10(x_kcat(i)),log10(x_kcat(i))],[log10(y_low(i)),log10(y_high(i))],'Color',[8,81,156]/255,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_low(i)),log10(y_low(i))],'Color',[8,81,156]/255,'LineWidth',0.5);
    line([log10(x_kcat(i))-0.1,log10(x_kcat(i))+0.1],[log10(y_high(i)),log10(y_high(i))],'Color',[8,81,156]/255,'LineWidth',0.5);
end
box on;
xlim([-2 6]);
ylim([-2 6]);
xticks([-2 0 2 4 6]);
yticks([-2 0 2 4 6]);
title(['R^2' '= ' num2str(round(RHO^2,3)) '; p' ' = ' num2str(round(PVAL,2)) '; N' ' = ' num2str(length(x_kcat))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) (/s)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[700 200 120 120]);
set(gca,'position',[0.2 0.2 0.7 0.7]);

figure();
histogram(log2(kapp4.fc));
set(gca,'FontSize',6,'FontName','Helvetica');
xlim([-1 12]);
xticks([0 5 10]);
ylim([0 50]);
xlabel('log2(max/min)','FontSize',7,'FontName','Helvetica');
ylabel('Count','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 120 120]);
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
    scatter(x,y,15,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
    
end
xlim([0 0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Mean kapp/kmax','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[100 100 120 120]);
set(gca,'position',[0.2 0.2 0.7 0.7]);


list_rho = zeros(size(list_satur,1),1);
list_p = zeros(size(list_satur,1),1);
for i = 1:size(list_satur,1)
    [m,n] = corr(list_gr',list_satur(i,:)','Type','Pearson');
    list_rho(i) = m;
    list_p(i) = n;
end
figure();
histogram(list_rho);
set(gca,'FontSize',6,'FontName','Helvetica');
title('Correlation between kapp/kmax and growth rate','FontSize',7,'FontName','Helvetica');
xlabel('Pearson r','FontSize',7,'FontName','Helvetica');
ylabel('Count','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 400 120 120]);
set(gca,'position',[0.2 0.2 0.7 0.7]);

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
h = boxplot([satur_low satur_high],'Notch','on','Symbol','.','OutlierSize',6,'Widths',0.3,'Colors',[8,81,156]/255);
set(h,{'linew'},{0.5});
xlim([0.5 2.5]);
xticks([1 2]);
set(gca, 'XColor','k');
set(gca, 'YColor','k');
set(gca,'XTickLabel',{['slow(N = ' num2str(lenlow) ')'],['fast(N = ' num2str(lenhigh) ')']});
set(gca,'FontSize',6,'FontName','Helvetica');
set(gca,'XColor','k');
set(gca,'YColor','k');
ylabel('kapp/kmax','FontSize',7,'FontName','Helvetica','Color','k');
% [~,ttestp] = ttest2(satur_low,satur_high,'Vartype','unequal');
ranksump = ranksum(satur_low,satur_high);
title(['rank sum test p' ' < 1e' num2str(ceil(log10(ranksump)))],'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[200 600 100 120]);
set(gca,'position',[0.25 0.2 0.7 0.7]);

%% kmax correlates with flux
figure();
[RHO,PVAL] = corr(log10(kapp4.pFBA),log10(kapp4.max),'Type','Pearson');
scatter(log10(kapp4.pFBA),log10(kapp4.max),15,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
box on;
xlim([-6 2]);
ylim([-4 6]);
% xticks([-2 0 2 4 6]);
% yticks([-2 0 2 4 6]);
title(['R^2' '= ' num2str(round(RHO^2,3)) '; p' ' < 1e' num2str(ceil(log10(PVAL))) '; N' ' = ' num2str(length(kapp4.pFBA))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (/s)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[400 400 120 120]);
set(gca,'position',[0.2 0.2 0.7 0.7]);

%% kapp comparison using estimated and measure fluxes
[num,txt,~] = xlsread('FluxMeasurements.xlsx','Relative');
rxnid = txt(2:end,1);
relfluxes = num;
titleid = txt(1,3:end);
clear num txt;
[qs,~,~] = xlsread('FluxMeasurements.xlsx','Absolute qs');
measuredfluxB = zeros(length(rxnid),0);
measuredfluxC = zeros(length(rxnid),0);
for i = 1:length(titleid)
    fluxtmp = relfluxes(:,i)*qs(i)/100;
    if contains(titleid{i},'Batch')
        measuredfluxB = [measuredfluxB fluxtmp];
    else
        measuredfluxC = [measuredfluxC fluxtmp];
    end
end
measuredfluxBmax = max(measuredfluxB,[],2);
measuredfluxCmax = max(measuredfluxC,[],2);

rxnlist_cmpr = cell(0,1);
kmax_estm = zeros(0,1);
kmax_meas = zeros(0,1);
for i = 1:length(rxnid)
    if ismember(rxnid(i),kapp4.rxn)
        if ismember(kapp4.max_cond(i),{'DiBartolomeo_GlucR1' 'DiBartolomeo_GlucR2' 'DiBartolomeo_GlucR3'})
            flux_estm = kapp.fluxes(ismember(kapp.rxn,rxnid(i)),ismember(kapp.condition,kapp4.max_cond(i)));
            flux_meas = measuredfluxBmax(i);
            kmax_estm = [kmax_estm;kapp4.max(i)];
            kmax_meas = [kmax_meas;kapp4.max(i)*flux_meas/flux_estm];
            rxnlist_cmpr = [rxnlist_cmpr;rxnid(i)];
        elseif ismember(kapp4.max_cond(i),{'Lahtvee_REF' 'Fluxes_Yu2_std_010'})
            flux_estm = kapp.fluxes(ismember(kapp.rxn,rxnid(i)),ismember(kapp.condition,kapp4.max_cond(i)));
            flux_meas = measuredfluxCmax(i);
            kmax_estm = [kmax_estm;kapp4.max(i)];
            kmax_meas = [kmax_meas;kapp4.max(i)*flux_meas/flux_estm];
            rxnlist_cmpr = [rxnlist_cmpr;rxnid(i)];
        end
    end
end
figure();
line([-2 6],[-2 6],'Color','k');
hold on;
[RHO,PVAL] = corr(log10(kmax_estm),log10(kmax_meas),'Type','Pearson');
scatter(log10(kmax_estm),log10(kmax_meas),15,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
box on;
xlim([-2 6]);
ylim([-2 6]);
xticks([-2 0 2 4 6]);
yticks([-2 0 2 4 6]);
set(gca,'FontSize',6,'FontName','Helvetica');
title('log10 kmax (in vivo) (/s)','FontSize',7,'FontName','Helvetica');
xlabel('using simulated flux','FontSize',7,'FontName','Helvetica');
ylabel('using measured flux','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[500 200 120 120]);
set(gca,'position',[0.2 0.2 0.7 0.7]);









