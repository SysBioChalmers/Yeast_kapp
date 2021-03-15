
load('ecYeastGEM_batch.mat'); %https://github.com/SysBioChalmers/ecModels
sol_org = optimizeCbModel(ecModel_batch);
maxmu_org = sol_org.f;
clear sol_org;

mu_list = 0.01:0.01:0.37;

model_org = changeObjective(ecModel_batch,'r_1714_REV');
model_org = changeRxnBounds(model_org,'r_1634',0,'b');% acetate production
model_org = changeRxnBounds(model_org,'r_2033',0,'b');% pyruvate production
model_org = changeRxnBounds(model_org,'r_1631',0,'b');% acetaldehyde production
model_org = changeRxnBounds(model_org,'r_1549',0,'b');%(R,R)-2,3-butanediol exchange

fluxes_org = zeros(length(model_org.rxns),length(mu_list));
for i = 1:length(mu_list)
    model_tmp = model_org;
    model_tmp = changeRxnBounds(model_tmp,'r_2111',mu_list(i),'b');
    sol = optimizeCbModel(model_tmp,'min');
    fluxes_org(:,i) = sol.x;
end

mu_org = fluxes_org(strcmp(model_org.rxns,'r_2111'),:);
glc_org = fluxes_org(strcmp(model_org.rxns,'r_1714_REV'),:);
etoh_org = fluxes_org(strcmp(model_org.rxns,'r_1761'),:);
o2_org = fluxes_org(strcmp(model_org.rxns,'r_1992_REV'),:);


%% replace kcat by kmax
load('kmax.mat');

model = ecModel_batch;

for i = 1:length(kapp3.rxn)
    rxntmp = kapp3.rxn{i};
    kmax = kapp3.max(i)*3600; % /h
    if length(rxntmp) == 6
        idxtmp = find(contains(model.rxns,rxntmp));
    elseif length(rxntmp) == 10
        if strcmp(rxntmp(8:end),'rvs')
            idxtmp = find(contains(model.rxns,rxntmp(1:6)) & contains(model.rxns,'REV'));
        elseif strcmp(rxntmp(8:end),'fwd')
            idxtmp = find(contains(model.rxns,rxntmp(1:6)) & ~contains(model.rxns,'REV') & ~contains(model.rxns,'arm'));
        end
    end
    for j = 1:length(idxtmp)
        metlist = model.mets(full(model.S(:,idxtmp(j))) < 0);
        protmet = metlist(contains(metlist,'prot_'));
        model.S(ismember(model.mets,protmet),idxtmp(j)) = -1/kmax;
    end
end

%% set model
sol_new = optimizeCbModel(model);
maxmu_new = sol_new.f;
clear sol_new;
prot_ub = model.ub(ismember(model.rxnNames,'prot_pool_exchange'));
model = changeRxnBounds(model,'prot_pool_exchange',prot_ub*maxmu_org/maxmu_new,'u');% decrease global saturation

model = changeObjective(model,'r_1714_REV');
model = changeRxnBounds(model,'r_1634',0,'b');% acetate production
model = changeRxnBounds(model,'r_2033',0,'b');% pyruvate production
model = changeRxnBounds(model,'r_1631',0,'b');% acetaldehyde production
model = changeRxnBounds(model,'r_1549',0,'b');%(R,R)-2,3-butanediol exchange

fluxes = zeros(length(model.rxns),length(mu_list));
for i = 1:length(mu_list)
    model_tmp = model;
    model_tmp = changeRxnBounds(model_tmp,'r_2111',mu_list(i),'b');
    sol = optimizeCbModel(model_tmp,'min');
    fluxes(:,i) = sol.x;
end

mu = fluxes(strcmp(model.rxns,'r_2111'),:);
glc = fluxes(strcmp(model.rxns,'r_1714_REV'),:);
etoh = fluxes(strcmp(model.rxns,'r_1761'),:);
o2 = fluxes(strcmp(model.rxns,'r_1992_REV'),:);

%% plot
% Yeast (PMID: 9603825)
fluxes_exp_yeast = [0.1  0.15  0.2  0.25  0.28  0.3   0.32  0.35  0.36  0.38 ; % mu
                    1.1  1.67  2.15 2.83  3.24  3.7   5.44  8.09  8.33  10.23; % glucose
                    0    0     0    0     0     0.51  4.42  6.91  6.71  14.91; % ethanol
                    2.73 2.5   5.07 6.8   8.3   8.8   6.83  6.6   7.1   4.19];% o2
figure();
subplot(1,2,1);
hold on;
box on;
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(2,:),'o','LineWidth',0.75,'Color',[55,126,184]/255,'MarkerSize',3);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(3,:),'o','LineWidth',0.75,'Color',[255,127,0]/255,'MarkerSize',3);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(4,:),'o','LineWidth',0.75,'Color',[77,175,74]/255,'MarkerSize',3);
plot(mu_org,glc_org,'-','LineWidth',0.75,'Color',[55,126,184]/255);
plot(mu_org,etoh_org,'-','LineWidth',0.75,'Color',[255,127,0]/255);
plot(mu_org,o2_org,'-','LineWidth',0.75,'Color',[77,175,74]/255);
xlim([0 0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
title('ecYeast with kcat','FontSize',7,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica');
legend({'Glucose uptake',...
        'Ethanol production'...
        'O2 uptake',},'FontSize',6,'FontName','Helvetica','location','nw');

subplot(1,2,2);
hold on;
box on;
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(2,:),'o','LineWidth',0.75,'Color',[55,126,184]/255,'MarkerSize',3);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(3,:),'o','LineWidth',0.75,'Color',[255,127,0]/255,'MarkerSize',3);
plot(fluxes_exp_yeast(1,:),fluxes_exp_yeast(4,:),'o','LineWidth',0.75,'Color',[77,175,74]/255,'MarkerSize',3);
plot(mu,glc,'-','LineWidth',0.75,'Color',[55,126,184]/255);
plot(mu,etoh,'-','LineWidth',0.75,'Color',[255,127,0]/255);
plot(mu,o2,'-','LineWidth',0.75,'Color',[77,175,74]/255);
xlim([0 0.4]);
set(gca,'FontSize',6,'FontName','Helvetica');
title('ecYeast with kmax','FontSize',7,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',7,'FontName','Helvetica');
ylabel('Flux (mmol/gCDW/h)','FontSize',7,'FontName','Helvetica');
legend({'Glucose uptake',...
        'Ethanol production'...
        'O2 uptake',},'FontSize',6,'FontName','Helvetica','location','nw');
set(gcf,'position',[100 400 240 100]);






