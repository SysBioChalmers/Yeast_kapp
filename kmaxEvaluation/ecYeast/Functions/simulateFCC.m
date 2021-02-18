function [enzyme_list,fcc_list] = simulateFCC(model)
% max growth 

model_ref = model;
sol_ref = optimizeCbModel(model_ref,'max');
mu_ref = sol_ref.f;

enzyme_list = model.enzNames;
fcc_list = zeros(length(enzyme_list),1);

for i = 1:length(model.enzymes)
    display([num2str(i),'/',num2str(length(model.enzymes))]);
    model_tmp = model;
    metidx = ismember(model_tmp.mets,strcat('prot_',model_tmp.enzymes{i}));
    rxnidx = full(model_tmp.S(metidx,:)) < 0;
    model_tmp.S(metidx,rxnidx) = model_tmp.S(metidx,rxnidx)/1.001;
    sol_tmp = optimizeCbModel(model_tmp,'max');
    mu_tmp = sol_tmp.f;
    fcc_list(i,1) = ((mu_tmp-mu_ref)/mu_ref)/0.001;
end

[fcc_list,idx] = sort(fcc_list,'ascend');
enzyme_list = enzyme_list(idx);


