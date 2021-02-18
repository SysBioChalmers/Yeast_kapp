function [sol_new, deviation, model_new] = searchFeasibleSol(model,exRxnList,exFluxes,osenseStr,step)
% This is for the case that no solution can be obtained for a given model
% with experimentally measured data as constraint. In order to get a
% feasible solution, we adjust exchange reaction rates gradually.

% first find optimization with the default objective
% second fix optimized objective and then minimize total enzymatic flux

model_org = changeRxnBounds(model,exRxnList,exFluxes,'b');
sol = optimizeCbModel(model_org,osenseStr);
if strcmp(sol.origStat,'OPTIMAL')
    model_adj = model_org;
    idxobj = model_org.c == 1;
    if sol.f >= 0 && strcmp(model_adj.osenseStr,'max')
        model_adj.ub(idxobj) = sol.f * 0.99999;
        model_adj.lb(idxobj) = sol.f * 0.99999;
    elseif sol.f >= 0 && strcmp(model_adj.osenseStr,'min')
        model_adj.ub(idxobj) = sol.f * 1.00001;
        model_adj.lb(idxobj) = sol.f * 1.00001;
    elseif sol.f < 0 && strcmp(model_adj.osenseStr,'max')
        model_adj.ub(idxobj) = sol.f * 1.00001;
        model_adj.lb(idxobj) = sol.f * 1.00001;
    elseif sol.f < 0 && strcmp(model_adj.osenseStr,'min')
        model_adj.ub(idxobj) = sol.f * 0.99999;
        model_adj.lb(idxobj) = sol.f * 0.99999;
    end
    [sol_new, model_new] = minimizeTotalFlux(model_adj);
    deviation = 0;
else
    deviation = step;
    while ~strcmp(sol.origStat,'OPTIMAL')
        exFluxes_lb = exFluxes;
        exFluxes_ub = exFluxes;

        exFluxes_lb(exFluxes_lb<0) = exFluxes_lb(exFluxes_lb<0)*(1+deviation);
        exFluxes_lb(exFluxes_lb>0) = exFluxes_lb(exFluxes_lb>0)*(1-deviation);

        exFluxes_ub(exFluxes_ub<0) = exFluxes_ub(exFluxes_ub<0)*(1-deviation);
        exFluxes_ub(exFluxes_ub>0) = exFluxes_ub(exFluxes_ub>0)*(1+deviation);

        model_tmp = model;
        model_tmp = changeRxnBounds(model_tmp,exRxnList,exFluxes_lb,'l');
        model_tmp = changeRxnBounds(model_tmp,exRxnList,exFluxes_ub,'u');
        sol = optimizeCbModel(model_tmp,osenseStr);

        deviation = deviation + step;
    end
    model_adj = model_tmp;
    idxobj = model_adj.c == 1;
    if sol.f >= 0 && strcmp(model_adj.osenseStr,'max')
        model_adj.ub(idxobj) = sol.f * 0.99999;
        model_adj.lb(idxobj) = sol.f * 0.99999;
    elseif sol.f >= 0 && strcmp(model_adj.osenseStr,'min')
        model_adj.ub(idxobj) = sol.f * 1.00001;
        model_adj.lb(idxobj) = sol.f * 1.00001;
    elseif sol.f < 0 && strcmp(model_adj.osenseStr,'max')
        model_adj.ub(idxobj) = sol.f * 1.00001;
        model_adj.lb(idxobj) = sol.f * 1.00001;
    elseif sol.f < 0 && strcmp(model_adj.osenseStr,'min')
        model_adj.ub(idxobj) = sol.f * 0.99999;
        model_adj.lb(idxobj) = sol.f * 0.99999;
    end
    [sol_new, model_new] = minimizeTotalFlux(model_adj);
end


