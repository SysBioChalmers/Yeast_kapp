%% blockRxns 
function model = removeTotFluxRxn(model,expID)

model.S = model.S(1:end-1,1:end-1);

model.rxns = model.rxns(1:end-1);
model.lb = model.lb(1:end-1);
model.ub = model.ub(1:end-1);
model.rxnGeneMat = model.rxnGeneMat(1:end-1,:);
model.c = model.c(1:end-1);
model.rules = model.rules(1:end-1);
model.rxnNames = model.rxnNames(1:end-1);
model.grRules = model.grRules(1:end-1);

model.mets = model.mets(1:end-1);
model.b = model.b(1:end-1);
model.csense = model.csense(1:end-1);
model.metFormulas = model.metFormulas(1:end-1);
model.metNames = model.metNames(1:end-1);


if contains(expID,{'Lahtvee','Yu2020_Clim','Yu2021_std_010','Yu2021_Gln_glc','Yu2021_Phe_std','Yu2021_Ile_std'})
    objrxn = 'r_4046'; % max ATP production
elseif contains(expID,'DiBartolomeo')
    objrxn = 'r_2111'; % max growth rate
elseif contains(expID,'Yu2021_Gln_N30')
    objrxn = 'r_1891'; % min gln uptake
elseif contains(expID,'Yu2021_Phe_N30')
    objrxn = 'r_1903'; % min phe uptake
elseif contains(expID,'Yu2021_Ile_N30')
    objrxn = 'r_1897'; % min ile uptake
else % the others are N-lim conditions
    objrxn = 'r_1654'; % min NH4 uptake
end

model = changeObjective(model, objrxn);
model.lb(ismember(model.rxns,objrxn)) = model.lb(ismember(model.rxns,objrxn)) - abs(model.lb(ismember(model.rxns,objrxn)))*0.01;
model.ub(ismember(model.rxns,objrxn)) = model.ub(ismember(model.rxns,objrxn)) + abs(model.ub(ismember(model.rxns,objrxn)))*0.01;






