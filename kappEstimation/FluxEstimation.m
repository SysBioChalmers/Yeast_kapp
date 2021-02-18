% load('GEM-yeast.mat');
% model_split = splitRevRxns(model);
% save('GEM-yeast-split.mat','model_split');

load('GEM-yeast-split.mat');
model = model_split;

model = setMedia(model);
model = blockRxns(model);

[num,txt,~] = xlsread('ProteomicsFlux.xlsx','Flux');
exRxnList = txt(2:end,1);
exFluxes = num;
expList = txt(1,3:end);
clear num txt;

for i = 1:length(expList)
    display([num2str(i),'/',num2str(length(expList))]);
    
    expID = expList{i};
        
    exFluxes_tmp = exFluxes(:,i)';
    exRxnList_tmp = exRxnList';
    
    idxnan = isnan(exFluxes_tmp);
    exFluxes_tmp = exFluxes_tmp(~idxnan);
    exRxnList_tmp = exRxnList_tmp(~idxnan);
    
    model_tmp = model;
    
    % change protein content for N-limited chemostats
    if contains(expID,'Yu1_CN30')
        model_tmp = scaleBioMass(model_tmp,'protein',0.3665,'carbohydrate',false);
    elseif contains(expID,'Yu1_CN50') || contains(expID,'Yu1_CN115')
        model_tmp = scaleBioMass(model_tmp,'protein',0.2635,'carbohydrate',false);
    elseif contains(expID,'Yu2_std_010')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.053444 0.050594 0.045962]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.35967375 0.377942893 0.336456714]),'carbohydrate',false);
    elseif contains(expID,'Yu2_N30_010')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.022625 0.02601 0.027257]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.209333929 0.221513357 0.245872214]),'carbohydrate',false);
    elseif contains(expID,'Yu2_N30_013')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.026366 0.025831 0.025297]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.3037245 0.256909821 0.266425]),'carbohydrate',false);
    elseif contains(expID,'Yu2_N30_018')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.039014 0.038658 0.037945]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.289261429 0.282029893 0.288119607]),'carbohydrate',false);
    elseif contains(expID,'Yu2_N30_030')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.047031 0.045606 0.052019]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.422473929 0.421332107 0.3996375]),'carbohydrate',false);
    elseif contains(expID,'Yu2_N30_035')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.073219 0.070546 0.088717]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.484154674 0.489921983 0.508864255]),'carbohydrate',false);
    elseif contains(expID,'Yu2_Gln_glc')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.071259 0.058789 0.064133]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.497834143 0.583851357 0.59412775]),'carbohydrate',false);
    elseif contains(expID,'Yu2_Gln_N30')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.026366 0.022268 0.017815]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.205527857 0.352822821 0.361957393]),'carbohydrate',false);
    elseif contains(expID,'Yu2_Phe_std')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.019596 0.020309 0.018171]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.425138179 0.520289964 0.452161286]),'carbohydrate',false);
    elseif contains(expID,'Yu2_Phe_N30')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.023337 0.024941 0.024406]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.417906643 0.378704107 0.398876286]),'carbohydrate',false);
    elseif contains(expID,'Yu2_Ile_std')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.029216 0.034739 0.030641]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.320851821 0.418667857 0.331508821]),'carbohydrate',false);
    elseif contains(expID,'Yu2_Ile_N30')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.024762 0.021912 0.024762]),'carbohydrate',false);
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.439981857 0.448355214 0.42361575]),'carbohydrate',false);
    end
    
    % set objective function for each condition
    if contains(expID,{'Lahtvee','Yu1_Clim','Yu2_std_010','Yu2_Gln_glc','Yu2_Phe_std','Yu2_Ile_std'})
        model_tmp = changeRxnBounds(model_tmp, 'r_4046', 0, 'l');
        model_tmp = changeRxnBounds(model_tmp, 'r_4046', 1000, 'u');
        model_tmp = changeObjective(model_tmp, 'r_4046'); % max ATP production
    elseif contains(expID,'DiBartolomeo')
        model_tmp = changeObjective(model_tmp, 'r_2111'); % max growth rate
    elseif contains(expID,'Yu2_Gln_N30')
        model_tmp = changeObjective(model_tmp, 'r_1891'); % min gln uptake
    elseif contains(expID,'Yu2_Phe_N30')
        model_tmp = changeObjective(model_tmp, 'r_1903'); % min phe uptake
    elseif contains(expID,'Yu2_Ile_N30')
        model_tmp = changeObjective(model_tmp, 'r_1897'); % min ile uptake
    else % the others are N-lim conditions
        model_tmp = changeObjective(model_tmp, 'r_1654'); % min NH4 uptake
    end
    
    [sol_new, deviation, model_new] = searchFeasibleSol(model_tmp,exRxnList_tmp,exFluxes_tmp,'max',0.001);
    
    display(['deviation: ',num2str(deviation)]);
    
    Fluxes.deviation = deviation;
    Fluxes.model = model_new;
    Fluxes.pFBA = sol_new.fluxes;
    cd Fluxes/;
    save(['Fluxes_' expID '.mat'],'Fluxes');
    cd ../;
end
