function [sol, model_new] = minimizeTotalFlux(model)
% The objective is to minimize the total flux through all enzymatic
% reactions. Should add a pseudo metabolite to each enzymatic reaction and
% an exchange reaction for the pseudo metabolite.

if ismember('r_tot_flux',model.rxns) % check if the pseudo metabolite is in the model
    model_new = model;
else % if it is a regular model but should be splitted into forward and reverse reactions
    model_tmp = model;
    model_tmp.rxns = [model_tmp.rxns;'r_tot_flux'];
    model_tmp.rxnNames = [model_tmp.rxnNames;'total flux'];
    model_tmp.rxnGeneMat = [model_tmp.rxnGeneMat;zeros(1,length(model_tmp.genes))];
    model_tmp.grRules = [model_tmp.grRules;{''}];
    model_tmp.rules = [model_tmp.rules;{''}];
    model_tmp.mets = [model_tmp.mets;'pseudo_met'];
    model_tmp.metNames = [model_tmp.metNames;'pseudo metabolite'];
    model_tmp.metFormulas = [model_tmp.metFormulas;{''}];
    
    model_tmp.b = [model_tmp.b;0];
    model_tmp.c = [zeros(length(model.rxns),1);1];
    model_tmp.lb = [model_tmp.lb;0];
    model_tmp.ub = [model_tmp.ub;100000];
    
    model_tmp.csense = [model_tmp.csense;'E'];
    model_tmp.osenseStr = 'min';
    
    model_tmp.S = [model_tmp.S zeros(size(model_tmp.S,1),1);zeros(1,length(model_tmp.rxns))];
    for i = 1:length(model_tmp.rxns)
        if ~ismember(model_tmp.grRules(i),'')
            model_tmp.S(ismember(model_tmp.mets,'pseudo_met'),i) = -1;
        end
    end
    model_tmp.S(ismember(model_tmp.mets,'pseudo_met'),ismember(model_tmp.rxns,'r_tot_flux')) = 1;
    model_new = model_tmp;
end

sol_tmp = optimizeCbModel(model_new,'min');

sol.origStat = sol_tmp.origStat;
sol.obj = sol_tmp.f;
sol.fluxes = sol_tmp.x;

