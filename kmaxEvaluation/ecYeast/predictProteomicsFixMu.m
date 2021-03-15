load('ecYeastGEM_batch.mat'); %https://github.com/SysBioChalmers/ecModels
sol_org = optimizeCbModel(ecModel_batch);

proteinlist = ecModel_batch.genes;

% list of carbon sources
CSlist = {'glucose' 'maltose' 'trehalose' 'fructose' 'sucrose' 'glycerol' 'acetate' 'pyruvate' 'lactate' 'oleate' 'galactose' 'raffinose'};
CSexlist = {'r_1714_REV' 'r_1931_REV' 'r_1650_REV' 'r_1709_REV' 'r_2058_REV' 'r_1808_REV' 'r_1634_REV' 'r_2033_REV' 'r_1551_REV' 'r_2189_REV' 'r_1710_REV' 'r_4043_REV'};
maxmulist = [0.46	0.14	0.15	0.46	0.46	0.24	0.17	0.23	0.11	0.16	0.34	0.33];

%% color
% colorveryhigh = [124,81,161]/255;
% colorhigh = [158,154,200]/255;
% colormedium = [218,218,235]/255;
% colorlow = [1,1,1];

colorveryhigh = [240,59,32]/255;
colorhigh = [254,178,76]/255;
colormedium = [255,237,160]/255;
colorlow = [1,1,1];

%% Assume a global turnover rate for all enzymes in ecYeast
model1 = ecModel_batch;
% Median of yeast (Saccharomyces cerevisiae) kcats in BRENDA is 7.84 /s.
% mdnkcat = 7.84 * 3600; % /h
% Median of Central-CE (PMID: 21506553)
mdnkcat = 79 * 3600; % /h
for i = 1:length(model1.rxns)
    if ~contains(model1.rxns{i},'draw_prot_')
        metlist = model1.mets(full(model1.S(:,i)) < 0);
        if any(contains(metlist,'prot_'))
            protmet = metlist(contains(metlist,'prot_'));
            model1.S(ismember(model1.mets,protmet),i) = -1/mdnkcat;
        end
    end
end
fluxes1 = simulateCarbonSourceFixMu(model1,CSexlist,maxmulist);
levellist1 = extractProteinLevels(model1,fluxes1,proteinlist);

%% Default ecYeast model
model2 = ecModel_batch;
fluxes2 = simulateCarbonSourceFixMu(model2,CSexlist,maxmulist);
levellist2 = extractProteinLevels(model2,fluxes2,proteinlist);

% load('kcat.mat');
% model2 = ecModel_batch;
% for i = 1:length(kcat.rxn)
%     rxntmp = kcat.rxn{i};
%     k = kcat.value(i)*3600; % /h
%     if length(rxntmp) == 6
%         idxtmp = find(contains(model2.rxns,rxntmp));
%     elseif length(rxntmp) == 10
%         if strcmp(rxntmp(8:end),'rvs')
%             idxtmp = find(contains(model2.rxns,rxntmp(1:6)) & contains(model2.rxns,'REV'));
%         elseif strcmp(rxntmp(8:end),'fwd')
%             idxtmp = find(contains(model2.rxns,rxntmp(1:6)) & ~contains(model2.rxns,'REV') & ~contains(model2.rxns,'arm'));
%         end
%     end
%     for j = 1:length(idxtmp)
%         metlist = model2.mets(full(model2.S(:,idxtmp(j))) < 0);
%         protmet = metlist(contains(metlist,'prot_'));
%         model2.S(ismember(model2.mets,protmet),idxtmp(j)) = -1/k;
%     end
% end
% fluxes2 = simulateCarbonSourceFixMu(model2,CSexlist,maxmulist);
% levellist2 = extractProteinLevels(model2,fluxes2,proteinlist);


%% Replace kcat of ecYeast by kmax

load('kmax.mat');
model4 = ecModel_batch;

for i = 1:length(kapp4.rxn)
    rxntmp = kapp4.rxn{i};
    kmax = kapp4.max(i)*3600; % /h
    if length(rxntmp) == 6
        idxtmp = find(contains(model4.rxns,rxntmp));
    elseif length(rxntmp) == 10
        if strcmp(rxntmp(8:end),'rvs')
            idxtmp = find(contains(model4.rxns,rxntmp(1:6)) & contains(model4.rxns,'REV'));
        elseif strcmp(rxntmp(8:end),'fwd')
            idxtmp = find(contains(model4.rxns,rxntmp(1:6)) & ~contains(model4.rxns,'REV') & ~contains(model4.rxns,'arm'));
        end
    end
    for j = 1:length(idxtmp)
        metlist = model4.mets(full(model4.S(:,idxtmp(j))) < 0);
        protmet = metlist(contains(metlist,'prot_'));
        model4.S(ismember(model4.mets,protmet),idxtmp(j)) = -1/kmax;
    end
end
fluxes4 = simulateCarbonSourceFixMu(model4,CSexlist,maxmulist);
levellist4 = extractProteinLevels(model4,fluxes4,proteinlist);


%% Replace kcat of ecYeast by mu-dependent kmax

load('kapp.mat');

[num,txt,~] = xlsread('ProteomicsFlux.xlsx','Flux');
exRxnList = txt(2:end,1);
exFluxes = num;
expList = txt(1,3:end);
clear num txt;
expGR = exFluxes(ismember(exRxnList,'r_2111'),:);
conditions = strrep(expList,'R1','');
conditions = strrep(conditions,'R2','');
conditions = strrep(conditions,'R3','');
condGR = zeros(1,length(kapp.condition));
for i = 1:length(kapp.condition)
    condGR(1,i) = mean(expGR(contains(conditions,kapp.condition(i))));
end

levellistmu = zeros(length(proteinlist),0);
for k = 1:length(CSlist)
    modelmu = ecModel_batch;
    mucutoff = maxmulist(k);
    idxcond = ~(condGR > mucutoff);
    for i = 1:length(kapp.rxn)
        rxntmp = kapp.rxn{i};
        kmax = max(kapp.values(i,idxcond))*3600; % /h
        if kmax > 0
            if length(rxntmp) == 6
                idxtmp = find(contains(modelmu.rxns,rxntmp));
            elseif length(rxntmp) == 10
                if strcmp(rxntmp(8:end),'rvs')
                    idxtmp = find(contains(modelmu.rxns,rxntmp(1:6)) & contains(modelmu.rxns,'REV'));
                elseif strcmp(rxntmp(8:end),'fwd')
                    idxtmp = find(contains(modelmu.rxns,rxntmp(1:6)) & ~contains(modelmu.rxns,'REV') & ~contains(modelmu.rxns,'arm'));
                end
            end
            for j = 1:length(idxtmp)
                metlist = modelmu.mets(full(modelmu.S(:,idxtmp(j))) < 0);
                protmet = metlist(contains(metlist,'prot_'));
                modelmu.S(ismember(modelmu.mets,protmet),idxtmp(j)) = -1/kmax;
            end
        end
    end
    fluxesmu = simulateCarbonSourceFixMu(modelmu,CSexlist(k),maxmulist(k));
    levellistmu = [levellistmu extractProteinLevels(modelmu,fluxesmu,proteinlist)];
end




%% Compared with exp data
% import exp data
% import molecular weight
[num,txt,~] = xlsread('UniProt.xlsx');
prot = txt(2:end,1);
MW = num;
clear num txt;
% glucose as reference
[num,txt,~] = xlsread('ProteomicsFlux.xlsx','DiBartolomeo');
gluc_protlist = txt(2:end,1);
condid = txt(1,2:end);
rawdata = mean(num(:,contains(condid,'Gluc')),2);
clear num txt;
gluc_abs = zeros(length(gluc_protlist),1);
for i = 1:length(gluc_protlist)
    if any(ismember(prot,gluc_protlist(i)))
        gluc_abs(i) = rawdata(i)*1000/MW(ismember(prot,gluc_protlist(i)));
    end
end
idxtmp = gluc_abs~=0;
gluc_abs = gluc_abs(idxtmp);
gluc_protlist = gluc_protlist(idxtmp);
clear idxtmp;

% [num,txt,~] = xlsread('ProteomicsCarbonSources.xlsx','RichAer');
% gluc_protlist = txt(2:end,1);
% rawdata = num(:,2);
% gluc_abs = zeros(length(gluc_protlist),1);
% for i = 1:length(gluc_protlist)
%     if any(ismember(prot,gluc_protlist(i)))
%         gluc_abs(i) = rawdata(i)*1000/MW(ismember(prot,gluc_protlist(i)));
%     end
% end
% clear num txt;

[num,txt,~] = xlsread('ProteomicsCarbonSources.xlsx','Paulo2015');
Paulo2015.protid = txt(2:end,1);
Paulo2015.csid = txt(1,2:end);
Paulo2015.data = num;
clear num txt;
[num,txt,~] = xlsread('ProteomicsCarbonSources.xlsx','Paulo2016');
Paulo2016.protid = txt(2:end,1);
Paulo2016.csid = txt(1,2:end);
Paulo2016.data = num;
clear num txt;

% analyze simulations
rmseCSlist = CSlist(~ismember(CSlist,'glucose'));
rmsedata = zeros(4,length(rmseCSlist));
ndata = cell(1,length(rmseCSlist));
for i = 1:length(rmseCSlist)
    % exp
    if ismember(rmseCSlist(i),{'galactose','raffinose'})
        protlisttmp = Paulo2015.protid;
        rawtmp = Paulo2015.data(:,ismember(Paulo2015.csid,rmseCSlist(i)));
        gluctmp = Paulo2015.data(:,ismember(Paulo2015.csid,'glucose'));
    else
        protlisttmp = Paulo2016.protid;
        rawtmp = Paulo2016.data(:,ismember(Paulo2016.csid,rmseCSlist(i)));
        gluctmp = Paulo2016.data(:,ismember(Paulo2016.csid,'glucose'));
    end
    
    fctmp = rawtmp./gluctmp;
    protlistexp = intersect(protlisttmp,gluc_protlist);
    [~,x] = ismember(protlistexp,gluc_protlist);
    [~,y] = ismember(protlistexp,protlisttmp);
    abslistexp = gluc_abs(x).*fctmp(y);
    idxtmp = abslistexp~=0;
    abslistexp = abslistexp(idxtmp);
    protlistexp = protlistexp(idxtmp);
    % sim
    idxtmptmp = ismember(CSlist,rmseCSlist(i));
    wholedata = [levellist1(:,idxtmptmp) levellist2(:,idxtmptmp) levellist4(:,idxtmptmp) levellistmu(:,idxtmptmp)];
    idxtmptmptmp = ~any(wholedata <= 0,2);
    protlistsim = proteinlist(idxtmptmptmp);
    abslistsim =  wholedata(idxtmptmptmp,:);
    % compare
    overlap = intersect(protlistexp,protlistsim);
    [~,m] = ismember(overlap,protlistexp);
    [~,n] = ismember(overlap,protlistsim);
    dataexp = abslistexp(m);
    datasim = abslistsim(n,:);
    
    dataexp = log10(dataexp);
    datasim = log10(datasim);
    
    rmsedata(:,i) = [sqrt(immse(datasim(:,1),dataexp))
                     sqrt(immse(datasim(:,2),dataexp))
                     sqrt(immse(datasim(:,3),dataexp))
                     sqrt(immse(datasim(:,4),dataexp))];
    ndata(1,i) = {['(N = ' num2str(length(dataexp)) ')']};
end

figure();
b = bar(1:length(rmseCSlist),rmsedata');
b(1).LineWidth = 0.5;
b(1).FaceColor = colorlow;
b(2).LineWidth = 0.5;
b(2).FaceColor = colormedium;
b(3).LineWidth = 0.5;
b(3).FaceColor = colorhigh;
b(4).LineWidth = 0.5;
b(4).FaceColor = colorveryhigh;
set(gca,'XTick',1:1:length(rmseCSlist));
set(gca,'XTickLabel',strcat(rmseCSlist,ndata));
ylim([0 2.5]);
legend({'global kcat' 'default kcat' 'kmax' 'kmax(mu)'},'Location','northwest','Orientation','horizontal','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('RMSE','FontSize',7,'FontName','Helvetica','Color','k');
title('Predictions of protein levels on various carbon sources','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[200 100 400 120]);
set(gca,'position',[0.1 0.2 0.85 0.7]);
box off;


