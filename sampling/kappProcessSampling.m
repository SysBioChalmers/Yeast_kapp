load('kapp_raw_sampling_median.mat');

conditions = strrep(kapp_raw_sampling.condition,'R1','');
conditions = strrep(conditions,'R2','');
conditions = strrep(conditions,'R3','');
conditions = unique(conditions);

kapp = struct();
kapp.condition = conditions;
kapp.rxn = kapp_raw_sampling.rxn;
kapp.protein = kapp_raw_sampling.protein;
kapp.values = zeros(length(kapp_raw_sampling.rxn),length(conditions));

for i = 1:length(conditions)
    idx = contains(kapp_raw_sampling.condition,conditions(i));
    values = kapp_raw_sampling.values(:,idx);
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

% add EC numbers
[~, rawUniprot, ~] = xlsread('ec_number.xlsx','Uniprot');
[~, rawKEGG, ~] = xlsread('ec_number.xlsx','KEGG');
id_Uniprot = rawUniprot(:,1);
ec_Uniprot = rawUniprot(:,2);
id_KEGG = rawKEGG(:,1);
ec_KEGG = rawKEGG(:,2);
id_list = unique([id_Uniprot;id_KEGG]);
ecdata = struct();
ecdata.id = cell(0,1);
ecdata.ec = cell(0,1);
for i = 1:length(id_list)
    id = id_list(i);
    if ismember(id,id_Uniprot)
        ec_U_tmp = ec_Uniprot(ismember(id_Uniprot,id));
        ec_U_tmp = split(ec_U_tmp);
        id_U_tmp = repelem(id,length(ec_U_tmp))';
        ecdata.id = [ecdata.id;id_U_tmp];
        ecdata.ec = [ecdata.ec;ec_U_tmp];
    end
    if ismember(id,id_KEGG)
        ec_K_tmp = ec_KEGG(ismember(id_KEGG,id));
        id_K_tmp = repelem(id,length(ec_K_tmp))';
        ecdata.id = [ecdata.id;id_K_tmp];
        ecdata.ec = [ecdata.ec;ec_K_tmp];
    end
end

kapp.EC = cell(length(kapp.rxn),1);
for i = 1:length(kapp.rxn)
    proteins = kapp.protein(i);
    proteins = strrep(proteins,'( ','');
    proteins = strrep(proteins,' )','');
    proteins = split(proteins,' or ');
    
    if any(ismember(proteins,ecdata.id))
        ectmp = cell(0,1);
        for j = 1:length(proteins)
            ectmp = [ectmp;ecdata.ec(ismember(ecdata.id,proteins(j)))];
        end
        ectmp = unique(ectmp);
        kapp.EC(i) = join(ectmp);
    else
        kapp.EC(i) = {''};
    end
end

[num,txt,~] = xlsread('ProteomicsFlux.xlsx','Flux');
exRxnList = txt(2:end,1);
exFluxes = num;
expList = txt(1,3:end);
clear num txt;
expGR = exFluxes(ismember(exRxnList,'r_2111'),:);
conditions = strrep(expList,'R1','');
conditions = strrep(conditions,'R2','');
conditions = strrep(conditions,'R3','');
kapp.condGR = zeros(1,length(kapp.condition));
for i = 1:length(kapp.condition)
    kapp.condGR(1,i) = mean(expGR(contains(conditions,kapp.condition(i))));
end

save('kapp_sampling_median.mat','kapp');
clear;

