load('kapp.mat');
load('kapp_raw.mat');

kapp4.max = zeros(0,1);
kapp4.rxn = cell(0,1);
kapp4.protein = cell(0,1);
kapp4.EC = cell(0,1);
kapp4.max_cond = cell(0,1);
kapp4.fc = zeros(0,1);
for i = 1:length(kapp.rxn)
    values_tmp = kapp.values(i,:);
    if sum(values_tmp > 0) >= 4 % && ~contains(kapp.protein(i),'or')
        kapp4.rxn = [kapp4.rxn;kapp.rxn(i)];
        kapp4.protein = [kapp4.protein;kapp.protein(i)];
        kapp4.EC = [kapp4.EC;kapp.EC(i)];
        kapp4.max = [kapp4.max;max(values_tmp)];
        values_raw = kapp_raw.values(ismember(kapp_raw.rxn,kapp.rxn(i)),:);
        kapp4.max_cond = [kapp4.max_cond;kapp_raw.condition(values_raw == max(values_tmp))];
        kapp4.fc = [kapp4.fc;max(values_tmp)/min(values_tmp)];
    end
end

kapp4.pFBA = zeros(0,1);
kapp4.maxFlux = zeros(0,1);
kapp4.minFlux = zeros(0,1);
kapp4.maxkmax = zeros(0,1);
kapp4.minkmax = zeros(0,1);

for i = 1:length(kapp4.rxn)
    display([num2str(i) '/' num2str(length(kapp4.rxn))]);
    rxnid = kapp4.rxn{i};
    load(['Fluxes_' kapp4.max_cond{i} '.mat']);
    pFBA_tmp = Fluxes.pFBA(ismember(Fluxes.model.rxns,rxnid));
    tot_flux_tmp = Fluxes.pFBA(ismember(Fluxes.model.rxns,'r_tot_flux'));
    
    model_tmp = Fluxes.model;
    model_tmp = changeRxnBounds(model_tmp,'r_tot_flux',tot_flux_tmp*0.9999,'l');
    model_tmp = changeRxnBounds(model_tmp,'r_tot_flux',tot_flux_tmp*1.0001,'u');
    model_tmp = changeObjective(model_tmp,rxnid);
    solmax = optimizeCbModel(model_tmp,'max');
    solmin = optimizeCbModel(model_tmp,'min');
    
    kapp4.pFBA = [kapp4.pFBA;round(pFBA_tmp,5)];
    kapp4.maxFlux = [kapp4.maxFlux;round(solmax.f,5)];
    kapp4.minFlux = [kapp4.minFlux;round(solmin.f,5)];
    kapp4.maxkmax = [kapp4.maxkmax;round(kapp4.max(i)*(round(solmax.f,5)/round(pFBA_tmp,5)),5)];
    kapp4.minkmax = [kapp4.minkmax;round(kapp4.max(i)*(round(solmin.f,5)/round(pFBA_tmp,5)),5)];
    
end

save('kmax.mat','kapp4');



