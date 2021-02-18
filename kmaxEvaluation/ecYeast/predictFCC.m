load('ecYeastGEM_batch.mat'); %https://github.com/SysBioChalmers/ecModels
sol_org = optimizeCbModel(ecModel_batch);
maxmu_org = sol_org.f;
clear sol_org;
prot_ub_org = ecModel_batch.ub(ismember(ecModel_batch.rxnNames,'prot_pool_exchange'));

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
[enzyme_list1,fcc_list1] = simulateFCC(model1);

%% Default ecYeast model
model2 = ecModel_batch;
[enzyme_list2,fcc_list2] = simulateFCC(model2);

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
[enzyme_list3,fcc_list3] = simulateFCC(model4);

%% plot
figure();
subplot(1,3,1);
x = fcc_list1;
y = enzyme_list1;
h = barh(1:length(y(end-19:end)),x(end-19:end),0.4,'FaceColor',[158,202,225]/255,'EdgeColor',[158,202,225]/255,'LineWidth',0.5);
xlim([0 0.12]);
set(gca,'YTick',1:1:length(y(end-19:end)));
set(gca,'YTickLabel',y(end-19:end));
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Flux control coefficient','FontSize',7,'FontName','Helvetica','Color','k');
title('ecYeast with a global kcat','FontSize',7,'FontName','Helvetica','Color','k');
subplot(1,3,2);
x = fcc_list2;
y = enzyme_list2;
h = barh(1:length(y(end-19:end)),x(end-19:end),0.4,'FaceColor',[66,146,198]/255,'EdgeColor',[66,146,198]/255,'LineWidth',0.5);
xlim([0 0.12]);
set(gca,'YTick',1:1:length(y(end-19:end)));
set(gca,'YTickLabel',y(end-19:end));
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Flux control coefficient','FontSize',7,'FontName','Helvetica','Color','k');
title('ecYeast with default kcat values','FontSize',7,'FontName','Helvetica','Color','k');
subplot(1,3,3);
x = fcc_list3;
y = enzyme_list3;
h = barh(1:length(y(end-19:end)),x(end-19:end),0.4,'FaceColor',[8,81,156]/255,'EdgeColor',[8,81,156]/255,'LineWidth',0.5);
xlim([0 0.12]);
set(gca,'YTick',1:1:length(y(end-19:end)));
set(gca,'YTickLabel',y(end-19:end));
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Flux control coefficient','FontSize',7,'FontName','Helvetica','Color','k');
title('ecYeast with kmax values','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[450 450 400 250]);

figure();
% [~,I1] = sort(enzyme_list1);
% fcclist_1 = fcc_list1(I1);
% [~,I2] = sort(enzyme_list2);
% fcclist_2 = fcc_list2(I2);
% [~,I3] = sort(enzyme_list3);
% fcclist_3 = fcc_list3(I3);
% 
% scatter(fcclist_1,fcclist_3)
% scatter(fcclist_1,fcclist_2)
% scatter(fcclist_2,fcclist_3)
% 
% [R,P] = corr(fcclist_2,fcclist_3,'Type','Pearson');

kmaxproteinlist = cell(0,1);
for i = 1:length(kapp4.protein)
    prottmp = kapp4.protein(i);
    if contains(prottmp,'or')
        prottmp = strrep(prottmp,'(','');
        prottmp = strrep(prottmp,')','');
        prottmp = strtrim(prottmp);
        prottmp = split(prottmp,' or ');
    end
    kmaxproteinlist = [kmaxproteinlist;prottmp];
end

ovlp = intersect(kmaxproteinlist,ecModel_batch.enzGenes);
[~,n] = ismember(ovlp,ecModel_batch.enzGenes);
kmaxenzlist = ecModel_batch.enzNames(n);

[~,n1] = ismember(kmaxenzlist,enzyme_list1);
fcclist1 = fcc_list1(n1);
[~,n2] = ismember(kmaxenzlist,enzyme_list2);
fcclist2 = fcc_list1(n2);
[~,n3] = ismember(kmaxenzlist,enzyme_list3);
fcclist3 = fcc_list1(n3);

subplot(3,1,1);
box on;
scatter(fcclist3,fcclist1,20,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
set(gca,'FontSize',6,'FontName','Helvetica');
title('FCC values of enzymes','FontSize',7,'FontName','Helvetica');
xlabel('kmax','FontSize',7,'FontName','Helvetica');
ylabel('global kcat','FontSize',7,'FontName','Helvetica');
subplot(3,1,2);
box on;
scatter(fcclist3,fcclist2,20,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('kmax','FontSize',7,'FontName','Helvetica');
ylabel('default kcat','FontSize',7,'FontName','Helvetica');
subplot(3,1,3);
box on;
scatter(fcclist1,fcclist2,20,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[8,81,156]/255,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('global kcat','FontSize',7,'FontName','Helvetica');
ylabel('default kcat','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 100 250]);


