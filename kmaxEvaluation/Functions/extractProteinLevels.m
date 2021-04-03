function levellist = extractProteinLevels(model,fluxes,proteinlist)

levellist = zeros(length(proteinlist),size(fluxes,2));
for i = 1:length(proteinlist)
    prottmp = proteinlist(i);
    if any(ismember(model.grRules,prottmp) & contains(model.rxns,'draw_prot_'))
        idxtmp = ismember(model.grRules,prottmp) & contains(model.rxns,'draw_prot_');
        levellist(i,:) = fluxes(idxtmp,:);
    end
end

