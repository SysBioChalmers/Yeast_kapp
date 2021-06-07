function FluxEstimationCluster(i)

% init cluster function %% please modify this part based on the position of
% those tool box
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'));
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'));
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'));
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/MATLAB-git'));
addpath(genpath('/cephyr/users/feiranl/Hebbe/kappEstimation'));

workdir = pwd;
cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox';
initCobraToolbox;
savepath '~/pathdef.m';
cd(workdir);


setRavenSolver('cobra');

load('GEM-yeast-split.mat');
model = model_split;

model = setMedia(model);
model = blockRxns(model);

[num,txt,~] = xlsread('ProteomicsFlux.xlsx','Flux');
exRxnList = txt(2:end,1);
exFluxes = num;
expList = txt(1,3:end);
clear num txt;


display([num2str(i),'/',num2str(length(expList))]);

expID = expList{i};

exFluxes_tmp = exFluxes(:,i)';
exRxnList_tmp = exRxnList';

idxnan = isnan(exFluxes_tmp);
exFluxes_tmp = exFluxes_tmp(~idxnan);
exRxnList_tmp = exRxnList_tmp(~idxnan);

model_tmp = model;

% change protein content for N-limited chemostats
if contains(expID,'Yu2020_CN30')
    model_tmp = scaleBioMass(model_tmp,'protein',0.3665,'carbohydrate',false);
elseif contains(expID,'Yu2020_CN50') || contains(expID,'Yu2020_CN115')
    model_tmp = scaleBioMass(model_tmp,'protein',0.2635,'carbohydrate',false);
elseif contains(expID,'Yu2021')
    [num,txt,~] = xlsread('ProteomicsFlux.xlsx','Yu2021RNAProtein');
    idtmp = txt(2:end,1);
    rna_cont = num(:,1);
    prot_cont = num(:,2);
    clear num txt;
    model_tmp = scaleBioMass(model_tmp,'RNA',rna_cont(ismember(idtmp,expID)),'carbohydrate',false);
    model_tmp = scaleBioMass(model_tmp,'protein',prot_cont(ismember(idtmp,expID)),'carbohydrate',false);
end

% set objective function for each condition
if contains(expID,{'Lahtvee','Yu2020_Clim','Yu2021_std_010','Yu2021_Gln_glc','Yu2021_Phe_std','Yu2021_Ile_std'})
    model_tmp = changeRxnBounds(model_tmp, 'r_4046', 0, 'l');
    model_tmp = changeRxnBounds(model_tmp, 'r_4046', 1000, 'u');
    model_tmp = changeObjective(model_tmp, 'r_4046'); % max ATP production
elseif contains(expID,'DiBartolomeo')
    model_tmp = changeObjective(model_tmp, 'r_2111'); % max growth rate
elseif contains(expID,'Yu2021_Gln_N30')
    model_tmp = changeObjective(model_tmp, 'r_1891'); % min gln uptake
elseif contains(expID,'Yu2021_Phe_N30')
    model_tmp = changeObjective(model_tmp, 'r_1903'); % min phe uptake
elseif contains(expID,'Yu2021_Ile_N30')
    model_tmp = changeObjective(model_tmp, 'r_1897'); % min ile uptake
else % the others are N-lim conditions
    model_tmp = changeObjective(model_tmp, 'r_1654'); % min NH4 uptake
end

[sol_new, deviation, model_new] = searchFeasibleSol(model_tmp,exRxnList_tmp,exFluxes_tmp,'max',0.001);

display(['deviation: ',num2str(deviation)]);

Fluxes.deviation = deviation;
Fluxes.model = model_new;
Fluxes.pFBA = sol_new.fluxes;

% random sampling
if Fluxes.deviation < 0.1
    model_sample = Fluxes.model;
    model_sample = removeTotFluxRxn(model_sample,expID);
    model_sample = ravenCobraWrapper(model_sample);
    [solutions,~] = randomSampling(model_sample,10000,1,0,1);
    Fluxes.SamplingFluxes = solutions;
%     tot_flux = Fluxes.pFBA(ismember(Fluxes.model.rxns,'r_tot_flux'));
%     model_sample = Fluxes.model;
%     model_sample = changeRxnBounds(model_sample,'r_tot_flux',tot_flux*0.9999,'l');
%     model_sample = changeRxnBounds(model_sample,'r_tot_flux',tot_flux*1.0001,'u');
%     model_sample = ravenCobraWrapper(model_sample);
%     [solutions,~] = randomSampling(model_sample,10000,0,0,1);
%     Fluxes.SamplingFluxes = solutions;
else
    Fluxes.SamplingFluxes = [];
end
cd Fluxes/;
save(['Fluxes_' expID '.mat'],'Fluxes');
cd ../;

end