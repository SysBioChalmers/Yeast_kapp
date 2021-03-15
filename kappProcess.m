load('kapp_raw.mat');

conditions = strrep(kapp_raw.condition,'R1','');
conditions = strrep(conditions,'R2','');
conditions = strrep(conditions,'R3','');
conditions = unique(conditions);

kapp = struct();
kapp.condition = conditions;
kapp.rxn = kapp_raw.rxn;
kapp.protein = kapp_raw.protein;
kapp.values = zeros(length(kapp_raw.rxn),length(conditions));

for i = 1:length(conditions)
    idx = contains(kapp_raw.condition,conditions(i));
    values = kapp_raw.values(:,idx);
    values(values == 0) = nan;
    if size(values,2) > 1
        for j = 1:size(values,1)
            cov = std(values(j,:),'omitnan')/mean(values(j,:),'omitnan');
            if cov > 0.5 || isnan(cov) % remove data with very high coefficient of variation
                kapp.values(j,i) = nan;
            else
                kapp.values(j,i) = max(values(j,:));
            end
        end
    elseif size(values,2) == 1
        kapp.values(:,i) = values;
    end
end
kapp.values(isnan(kapp.values)) = 0;

save('kapp.mat','kapp');
clear;

