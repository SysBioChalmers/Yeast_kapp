CalType = 'median';
% CalType = 'mean';


load(['kapp_sampling_' CalType '.mat']);
load(['kapp_raw_sampling_' CalType '.mat']);

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
        values_raw = kapp_raw_sampling.values(ismember(kapp_raw_sampling.rxn,kapp.rxn(i)),:);
        kapp4.max_cond = [kapp4.max_cond;kapp_raw_sampling.condition(values_raw == max(values_tmp))];
        kapp4.fc = [kapp4.fc;max(values_tmp)/min(values_tmp)];
    end
end

save(['kmax_sampling_' CalType '.mat'],'kapp4');



