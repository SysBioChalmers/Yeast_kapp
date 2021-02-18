load('ecYeastGEM_batch.mat'); %https://github.com/SysBioChalmers/ecModels
sol_org = optimizeCbModel(ecModel_batch);
maxmu_org = sol_org.f;
clear sol_org;
prot_ub_org = ecModel_batch.ub(ismember(ecModel_batch.rxnNames,'prot_pool_exchange'));

proteinlist = ecModel_batch.genes;

% list of carbon sources
CSlist = {'glucose' 'maltose' 'trehalose' 'fructose' 'sucrose' 'glycerol' 'acetate' 'pyruvate' 'lactate' 'oleate' 'galactose' 'raffinose'};
CSexlist = {'r_1714_REV' 'r_1931_REV' 'r_1650_REV' 'r_1709_REV' 'r_2058_REV' 'r_1808_REV' 'r_1634_REV' 'r_2033_REV' 'r_1551_REV' 'r_2189_REV' 'r_1710_REV' 'r_4043_REV'};

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
sol1 = optimizeCbModel(model1);
maxmu1 = sol1.f;
clear sol1;
model1 = changeRxnBounds(model1,'prot_pool_exchange',prot_ub_org*maxmu_org/maxmu1,'u');% adjust global saturation
% writeCbModel(model1,'xls','model1.xls');
fluxes1 = simulateCarbonSource(model1,CSexlist);
levellist1 = extractProteinLevels(model1,fluxes1,proteinlist);

%% Default ecYeast model
model2 = ecModel_batch;
fluxes2 = simulateCarbonSource(model2,CSexlist);
levellist2 = extractProteinLevels(model2,fluxes2,proteinlist);

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

sol3 = optimizeCbModel(model4);
maxmu3 = sol3.f;
clear sol3;
model4 = changeRxnBounds(model4,'prot_pool_exchange',prot_ub_org*maxmu_org/maxmu3,'u');% adjust global saturation
% writeCbModel(model4,'xls','model4.xls');
fluxes3 = simulateCarbonSource(model4,CSexlist);
levellist3 = extractProteinLevels(model4,fluxes3,proteinlist);

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
rmsedata = zeros(3,length(rmseCSlist));
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
    wholedata = [levellist1(:,idxtmptmp) levellist2(:,idxtmptmp) levellist3(:,idxtmptmp)];
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
                     sqrt(immse(datasim(:,3),dataexp))];
    ndata(1,i) = {['(N = ' num2str(length(dataexp)) ')']};
end

figure();
b = bar(1:length(rmseCSlist),rmsedata');
b(1).LineWidth = 0.1;
b(1).FaceColor = [158,202,225]/255;
b(2).LineWidth = 0.1;
b(2).FaceColor = [66,146,198]/255;
b(3).LineWidth = 0.1;
b(3).FaceColor = [8,81,156]/255;
set(gca,'XTick',1:1:length(rmseCSlist));
set(gca,'XTickLabel',strcat(rmseCSlist,ndata));
ylim([0 2.5]);
legend({'global kcat' 'default kcat' 'kmax'},'Location','northwest','Orientation','horizontal','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('RMSE','FontSize',7,'FontName','Helvetica','Color','k');
title('Predictions of protein levels on various carbon sources','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[200 100 360 120]);
set(gca,'position',[0.1 0.2 0.85 0.7]);
box off;

% mu_list = 0.01:0.01:0.37;
% model_test = model4;
% model_test = changeObjective(model_test,'r_1714_REV');
% model_test = changeRxnBounds(model_test,'r_1634',0,'b');% acetate production
% model_test = changeRxnBounds(model_test,'r_2033',0,'b');% pyruvate production
% model_test = changeRxnBounds(model_test,'r_1631',0,'b');% acetaldehyde production
% model_test = changeRxnBounds(model_test,'r_1549',0,'b');%(R,R)-2,3-butanediol exchange
% fluxes_org = zeros(length(model_test.rxns),length(mu_list));
% for i = 1:length(mu_list)
%     model_tmp = model_test;
%     model_tmp = changeRxnBounds(model_tmp,'r_2111',mu_list(i),'b');
%     sol = optimizeCbModel(model_tmp,'min');
%     fluxes_org(:,i) = sol.x;
% end
% mu_org = fluxes_org(strcmp(model_test.rxns,'r_2111'),:);
% glc_org = fluxes_org(strcmp(model_test.rxns,'r_1714_REV'),:);
% etoh_org = fluxes_org(strcmp(model_test.rxns,'r_1761'),:);
% o2_org = fluxes_org(strcmp(model_test.rxns,'r_1992_REV'),:);
% figure();
% hold on;
% box on;
% plot(mu_org,glc_org,'-','LineWidth',0.75,'Color',[55,126,184]/255);
% plot(mu_org,etoh_org,'-','LineWidth',0.75,'Color',[255,127,0]/255);
% plot(mu_org,o2_org,'-','LineWidth',0.75,'Color',[77,175,74]/255);



