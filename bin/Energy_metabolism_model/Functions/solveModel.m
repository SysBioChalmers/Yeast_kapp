function sol = solveModel(model,objective,osenseStr,prot_cost_info,tot_prot_weight,factor_hy,factor_ly)


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

% Determine objective
model.c(strcmp(model.rxns,objective)) = 1;

% Add protein cost infomation
[nMets,nRxns] = size(model.S);
cost_list = zeros(1,nRxns);
for i = 1:nRxns
    rxnid = model.rxns{i};
    if contains(rxnid,'_HY')
        id_tmp = rxnid(1:end-3);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp))/factor_hy;
        else
            cost = 0;
        end
    elseif contains(rxnid,'_LY')
        id_tmp = rxnid(1:end-3);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp))/factor_ly;
        else
            cost = 0;
        end
    elseif contains(rxnid,'_Bio')
        id_tmp = rxnid(1:end-4);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp));
        else
            cost = 0;
        end
    else
        id_tmp = rxnid;
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp));
        else
            cost = 0;
        end
    end
    cost_list(1,i) = cost;
end

% Construct LP
A = cost_list;
b = tot_prot_weight*1000;

Aeq = model.S;
beq = zeros(nMets,1);

lb = model.lb;
ub = model.ub;

%linprog always runs minimization
if osenseStr == 'max'
    f = -model.c;
elseif osenseStr == 'min'
    f = model.c;
end

[x,~,exitflag,~] = linprog(f,A,b,Aeq,beq,lb,ub);

sol = struct();
sol.fluxes = x;
sol.exitflag = exitflag;

if exitflag == 1
	sol.protUsage = cost_list * x / 1000;
end

    