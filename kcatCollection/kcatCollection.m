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
        if kcat.value(idxtmp) == kcattmp
            kcat.ref(idxtmp) = strcat(kcat.ref(idxtmp),{'; '},reftmp);
        elseif kcat.value(idxtmp) < kcattmp
            kcat.value(idxtmp) = kcattmp;
            kcat.ref(idxtmp) = reftmp;
            kcat.HeterExp(idxtmp) = hetmp;
        end
    end
end

cd ../;
save('kcat.mat','kcat');
clear;
