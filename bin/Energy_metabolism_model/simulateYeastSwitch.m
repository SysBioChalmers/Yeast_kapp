% Simulate yeast metabolic switch
% Search for the ratio of LY/HY for best-fit simulations

% Related to Fig 3BD
% Related to Fig S2B

%% Yeast CEN.PK 113-7D

% Chemostat data (DOI: 10.1038/s42255-018-0006-7)
% Data were obtained from Fig 3 of (DOI: 10.1038/s42255-018-0006-7), which 
% were originally reported in (PMID: 21354323) and (PMID: 9603825).
chm_data = [0.02	0.05	0.05	0.10	0.15	0.20	0.20	0.25	0.28	0.30	0.30	0.32	0.35	0.35	0.36	0.38	0.38    % Growth rate
            0.25	0.57	0.60	1.10	1.67	2.15	2.22	2.83	3.24	3.70	4.67	5.44	7.62	8.09	8.33	10.23	13.18   % Glucose rate
            0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.51	3.24	4.42	8.72	6.91	6.71	14.91	12.50]; % Ethanol rate

% Batch data
bch_data = [0.40        0.37        0.46  % Growth rate
            19.9        19.6        22.4  % Glucose rate
            29.6        30.5        36.1];% Ethanol rate
% PMID:     19684065    11157958    20199578
exp_data = [chm_data bch_data];

model_yeast = xls2model('Model_yeast.xlsx');

prot_cost_yeast = struct();
[num, txt, ~] = xlsread('Model_yeast.xlsx','Protein_cost_info');
prot_cost_yeast.id = txt(2:end,1);
prot_cost_yeast.value = num;
clear num txt;

mu_list = 0.01:0.005:0.5;

f_hy_yeast = 1;
f_ly_yeast = 1;

% Estimate protein allocation
prot_list = zeros(1,length(bch_data));
for i = 1:length(bch_data)
    model = model_yeast;
    model = changeRxnBounds(model, 'EXbiomass', bch_data(1,i), 'b');
    model = changeRxnBounds(model, 'EXglc', -bch_data(2,i), 'b');
    model = changeRxnBounds(model, 'EXetoh', bch_data(3,i), 'b');
    [minProt, fluxes] = minimizeProtein(model, prot_cost_yeast, 'ATPM', 'max', f_hy_yeast,f_ly_yeast);
    if ~isempty(fluxes)
        prot_list(1,i) = minProt;
    else
        [minProt, ~] = obtainFeasibleSol(model,prot_cost_yeast,'Yeast',0.01,'ATPM','max',f_hy_yeast,f_ly_yeast);
        prot_list(1,i) = minProt;
    end
end
min_prot_yeast = median(prot_list);
clear fluxes minProt model prot_list i;

% Generate figure

fluxes_sim_yeast = zeros(3,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);
    model = changeRxnBounds(model_yeast, 'EXbiomass', mu, 'b');
    sol = solveModel(model,'EXglc','max',prot_cost_yeast,min_prot_yeast,f_hy_yeast,f_ly_yeast);
    if sol.exitflag == 1
        glc = -sol.fluxes(strcmp(model.rxns,'EXglc'));
        eth = sol.fluxes(strcmp(model.rxns,'EXetoh'));
        fluxes_sim_yeast(:,i) = [mu; glc; eth];
    end
end
fluxes_sim_yeast = fluxes_sim_yeast(:,1:length(find(fluxes_sim_yeast(1,:))));
clear eth glc i min_prot_list mu sol;

figure('Name','yeast_1');
hold on;
box on;
plot(chm_data(2,:),chm_data(3,:),'o','LineWidth',0.75,'Color','k','MarkerSize',8);
plot(bch_data(2,:),bch_data(3,:),'^','LineWidth',0.75,'Color','k','MarkerSize',10);
plot(fluxes_sim_yeast(2,:),fluxes_sim_yeast(3,:),'-','LineWidth',0.75,'Color','k');

xlim([0 25]);

set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Glucose uptake (mmol/gCDW/h)','FontSize',14,'FontName','Helvetica');
ylabel(['Ethanol production',char(13,10)','(mmol/gCDW/h)'],'FontSize',14,'FontName','Helvetica');

legend({'Chemostat data',...
        'Batch data',...
        'Predicted data'},'FontSize',12,'FontName','Helvetica','location','nw');


set(gcf,'position',[300 400 240 185]);
set(gca,'position',[0.2 0.18 0.76 0.8]);

clear bch_data chm_data exp_data fluxes_sim_yeast model mu_list;
clear f_hy_yeast f_ly_yeast min_prot_yeast;
