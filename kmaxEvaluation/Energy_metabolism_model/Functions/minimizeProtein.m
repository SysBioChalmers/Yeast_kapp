function [minProt, fluxes] = minimizeProtein(model, prot_cost_info, objective, osenseStr, factor_hy, factor_ly, step_min)

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

prot_low = 0;
prot_high = 1;

while prot_high-prot_low > step_min
    prot_mid = (prot_low+prot_high)/2;
	
    sol = solveModel(model,objective,osenseStr,prot_cost_info,prot_mid,factor_hy,factor_ly);
    
	if sol.exitflag == 1
        prot_high = prot_mid;
    else
        prot_low = prot_mid;
	end
end

minProt = prot_high;

sol = solveModel(model,'EXbiomass','max',prot_cost_info,minProt);
fluxes = sol.fluxes;
