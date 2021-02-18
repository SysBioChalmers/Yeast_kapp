function [minProt, fluxes] = obtainFeasibleSol(model,prot_cost_info,orgName,step,objective,osenseStr,factor_hy,factor_ly,step_min)
% This is for the case that no solution can be obtained for a given model
% with experimentally measured data as constraint. In order to get a
% feasible solution, we adjust exchange reaction rates gradually.
if exist('objective', 'var')
    if isempty(objective)
        objective = 'EXbiomass';
    end
else
    objective = 'EXbiomass';
end
    
if exist('osenseStr', 'var')
    if isempty(osenseStr)
        osenseStr = 'max';
    end
else
    osenseStr = 'max';
end

if exist('factor_hy', 'var')
    if isempty(factor_hy)
        factor_hy = 1;
    end
else
    factor_hy = 1;
end

if exist('factor_ly', 'var')
    if isempty(factor_ly)
        factor_ly = 1;
    end
else
    factor_ly = 1;
end

if exist('step_min', 'var')
    if isempty(step_min)
        step_min = 0.0001;
    end
else
    step_min = 0.0001;
end

% bio = 'EXbiomass';
glc = 'EXglc';

if strcmp(orgName, 'Ecoli')
    byprod = 'EXac';
elseif strcmp(orgName, 'Yeast')
    byprod = 'EXetoh';
end

[minProt, fluxes] = minimizeProtein(model,prot_cost_info,objective,osenseStr,factor_hy,factor_ly,step_min);

while isempty(fluxes)
    
    glc_lb = model.lb(strcmp(model.rxns,glc));
    glc_ub = model.ub(strcmp(model.rxns,glc));
%     bio_lb = model.lb(strcmp(model.rxns,bio));
%     bio_ub = model.ub(strcmp(model.rxns,bio));
    byprod_lb = model.lb(strcmp(model.rxns,byprod));
    byprod_ub = model.ub(strcmp(model.rxns,byprod));
    
    model = changeRxnBounds(model, glc, glc_lb*(1+step), 'l');
    model = changeRxnBounds(model, glc, glc_ub*(1-step), 'u');
%     model = changeRxnBounds(model, bio, bio_lb*(1-step), 'l');
%     model = changeRxnBounds(model, bio, bio_ub*(1+step), 'u');
    model = changeRxnBounds(model, byprod, byprod_lb*(1-step), 'l');
    model = changeRxnBounds(model, byprod, byprod_ub*(1+step), 'u');
    
    [minProt, fluxes] = minimizeProtein(model,prot_cost_info,objective,osenseStr,factor_hy,factor_ly,step_min);
    
    step = step * 2;
end





