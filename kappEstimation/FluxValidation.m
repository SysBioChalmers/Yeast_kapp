%% simulated flux data
% glucose batch
load('Fluxes_DiBartolomeo_GlucR1.mat');
rxnlist = Fluxes.model.rxns;
fluxB1 = Fluxes.pFBA;
load('Fluxes_DiBartolomeo_GlucR2.mat');
fluxB2 = Fluxes.pFBA;
load('Fluxes_DiBartolomeo_GlucR3.mat');
fluxB3 = Fluxes.pFBA;
fluxB = mean([fluxB1 fluxB2 fluxB3],2);
relfluxB = fluxB*100/-fluxB(ismember(rxnlist,'r_1714'));

% glucose-limited chemostat D=0.1/h
load('Fluxes_Yu2_std_010R1.mat');
fluxC1 = Fluxes.pFBA;
load('Fluxes_Yu2_std_010R2.mat');
fluxC2 = Fluxes.pFBA;
load('Fluxes_Yu2_std_010R3.mat');
fluxC3 = Fluxes.pFBA;
fluxC = mean([fluxC1 fluxC2 fluxC3],2);
relfluxC = fluxC*100/-fluxC(ismember(rxnlist,'r_1714'));

%% measured flux data
[num,txt,~] = xlsread('FluxMeasurements.xlsx','Relative');
rxnid = txt(2:end,1);
fluxes = num;
titleid = txt(1,3:end);
clear num txt;

%% figure

% glucose batch
figure();
idx = contains(titleid,'Batch');
simulatedFlux = relfluxB;
titletmp = titleid(idx);
for i = 1:length(titletmp)
    measured = fluxes(:,ismember(titleid,titletmp{i}));
    idxnan = isnan(measured);
    rxnidtmp = rxnid(~idxnan);
    measured = measured(~idxnan);
    [~,b] = ismember(rxnidtmp,rxnlist);
    simulated = simulatedFlux(b);
    [RHOtmp,PVALtmp] = corr(simulated,measured,'Type','Pearson');
    subplot(1,length(titletmp),i);
    line([0 100],[0 100],'Color','k');
    hold on;
    box on;
    scatter(simulated,measured,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
    xlim([0 100]);
    ylim([0 100]);
    set(gca,'FontSize',6,'FontName','Helvetica');
    xlabel('Simulated flux','FontSize',7,'FontName','Helvetica');
    ylabel('Measured flux','FontSize',7,'FontName','Helvetica');
    titlename = strrep(titletmp{i},'_',' ');
    titlename = titlename(1:strfind(titlename,'(')-2);
	title(titlename,'FontSize',6,'FontName','Helvetica','FontWeight','bold');
    text(10,90,['R^2' '= ' num2str(round(RHOtmp^2,2))],'FontSize',6,'FontName','Helvetica');
    text(10,70,['p < 1e' num2str(ceil(log10(PVALtmp)))],'FontSize',6,'FontName','Helvetica');
end
set(gcf,'position',[400 400 100*length(titletmp) 80]);


% glucose-limited chemostat D=0.1/h
figure();
idx = contains(titleid,'ClimD01');
simulatedFlux = relfluxC;
titletmp = titleid(idx);
for i = 1:length(titletmp)
    measured = fluxes(:,ismember(titleid,titletmp{i}));
    idxnan = isnan(measured);
    rxnidtmp = rxnid(~idxnan);
    measured = measured(~idxnan);
    [~,b] = ismember(rxnidtmp,rxnlist);
    simulated = simulatedFlux(b);
    [RHOtmp,PVALtmp] = corr(simulated,measured,'Type','Pearson');
    subplot(1,length(titletmp),i);
    line([0 100],[0 100],'Color','k');
    hold on;
    box on;
    scatter(simulated,measured,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',1,'MarkerEdgeAlpha',0);
    xlim([0 100]);
    ylim([0 100]);
    set(gca,'FontSize',6,'FontName','Helvetica');
    xlabel('Simulated flux','FontSize',7,'FontName','Helvetica');
    ylabel('Measured flux','FontSize',7,'FontName','Helvetica');
    titlename = strrep(titletmp{i},'_',' ');
    titlename = titlename(1:strfind(titlename,'(')-2);
	title(titlename,'FontSize',6,'FontName','Helvetica','FontWeight','bold');
    text(10,90,['R^2' '= ' num2str(round(RHOtmp^2,2))],'FontSize',6,'FontName','Helvetica');
    text(10,70,['p < 1e' num2str(ceil(log10(PVALtmp)))],'FontSize',6,'FontName','Helvetica');
end
set(gcf,'position',[300 400 100*length(titletmp) 80]);







