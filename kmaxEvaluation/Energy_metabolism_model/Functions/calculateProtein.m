function [HY, LY, Bio] = calculateProtein(model, prot_cost_info, fluxes, factor_hy, factor_ly)

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


HY = 0;
LY = 0;
Bio = 0;

for i = 1:length(model.rxns)
    rxnid = model.rxns{i};
    if contains(rxnid,'_HY')
        id_tmp = rxnid(1:end-3);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp))/factor_hy;
            flux = fluxes(i);
            HY = HY + cost * flux;
        end
        
        
    elseif contains(rxnid,'_LY')
        id_tmp = rxnid(1:end-3);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp))/factor_ly;
            flux = fluxes(i);
            LY = LY + cost * flux;
        end
        
        
    elseif contains(rxnid,'_Bio')
        id_tmp = rxnid(1:end-4);
        if any(strcmp(prot_cost_info.id,id_tmp))
            cost = prot_cost_info.value(strcmp(prot_cost_info.id,id_tmp));
            flux = fluxes(i);
            Bio = Bio + cost * flux;
        end
    end
end

HY = HY / 1000;
LY = LY / 1000;
Bio = Bio / 1000;