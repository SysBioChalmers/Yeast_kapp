[num,txt,~] = xlsread('kcatCollected.xlsx','BRENDA20210203');
rxnB = txt(2:end,1);
refB = txt(2:end,3);
kcatB = num(:,1);
HeB = num(:,5);

[num,txt,~] = xlsread('kcatCollected.xlsx','Literature');
rxnL = txt(2:end,1);
refL = txt(2:end,3);
kcatL = num(:,1);
HeL = num(:,4);

[num,txt,~] = xlsread('kcatCollected.xlsx','SABIORK20210203');
rxnS = txt(2:end,1);
refS = txt(2:end,3);
kcatS = num(:,1);
HeS = num(:,8);

[num,txt,~] = xlsread('kcatCollected.xlsx','SA');
rxnSA = txt(2:end,1);
refSA = txt(2:end,3);
kcatSA = num(:,1);
HeSA = num(:,6);

allrxn = [rxnB;rxnL;rxnS;rxnSA];
allref = [refB;refL;refS;refSA];
allkcat = [kcatB;kcatL;kcatS;kcatSA];
allHeterExp = [HeB;HeL;HeS;HeSA];

kcat = struct();
kcat.value = zeros(0,1);
kcat.rxn = cell(0,1);
kcat.ref = cell(0,1);
kcat.HeterExp = zeros(0,1);

for i = 1:length(allrxn)
    rxntmp = allrxn(i);
    kcattmp = allkcat(i);
    reftmp = allref(i);
    hetmp = allHeterExp(i);
    if ~ismember(rxntmp,kcat.rxn)
        kcat.rxn = [kcat.rxn;rxntmp];
        kcat.value = [kcat.value;kcattmp];
        kcat.ref = [kcat.ref;reftmp];
        kcat.HeterExp = [kcat.HeterExp;hetmp];
    else
        idxtmp = ismember(kcat.rxn,rxntmp);
        if kcat.HeterExp(idxtmp) == 1 && hetmp ~= 1
            kcat.value(idxtmp) = kcattmp;
            kcat.ref(idxtmp) = reftmp;
            kcat.HeterExp(idxtmp) = hetmp;
        elseif (kcat.HeterExp(idxtmp) == 1 && hetmp == 1) || (kcat.HeterExp(idxtmp) ~= 1 && hetmp ~= 1)
            if kcat.value(idxtmp) == kcattmp
                kcat.ref(idxtmp) = strcat(kcat.ref(idxtmp),{'; '},reftmp);
            elseif kcat.value(idxtmp) < kcattmp
                kcat.value(idxtmp) = kcattmp;
                kcat.ref(idxtmp) = reftmp;
                kcat.HeterExp(idxtmp) = hetmp;
            end
        end
    end
end



% add proteins and EC numbers
load('GEM-yeast-split.mat');

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

kcat.EC = cell(length(kcat.rxn),1);
kcat.protein = cell(length(kcat.rxn),1);

for i = 1:length(kcat.rxn)
    
    proteins = model_split.grRules(ismember(model_split.rxns,kcat.rxn(i)));
        kcat.protein(i) = proteins;
        proteins = strrep(proteins,'( ','');
        proteins = strrep(proteins,' )','');
        proteins = split(proteins,' or ');

        if any(ismember(proteins,ecdata.id))
            ectmp = cell(0,1);
            for j = 1:length(proteins)
                ectmp = [ectmp;ecdata.ec(ismember(ecdata.id,proteins(j)))];
            end
            ectmp = unique(ectmp);
            kcat.EC(i) = join(ectmp);
        else
            kcat.EC(i) = {''};
        end
end


cd ../;
save('kcat.mat','kcat');
clear;
