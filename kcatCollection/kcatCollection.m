[num,txt,~] = xlsread('kcatCollected.xlsx','BRENDA20210203');
rxnB = txt(2:end,1);
refB = txt(2:end,3);
kcatB = num;
[num,txt,~] = xlsread('kcatCollected.xlsx','Literature');
rxnL = txt(2:end,1);
refL = txt(2:end,3);
kcatL = num;
[num,txt,~] = xlsread('kcatCollected.xlsx','SABIORK20210203');
rxnS = txt(2:end,1);
refS = txt(2:end,3);
kcatS = num;

allrxn = [rxnB;rxnL;rxnS];
allref = [refB;refL;refS];
allkcat = [kcatB;kcatL;kcatS];

kcat = struct();
kcat.value = zeros(0,1);
kcat.rxn = cell(0,1);
kcat.ref = cell(0,1);

for i = 1:length(allrxn)
    rxntmp = allrxn(i);
    kcattmp = allkcat(i);
    reftmp = allref(i);
    if ~ismember(rxntmp,kcat.rxn)
        kcat.rxn = [kcat.rxn;rxntmp];
        kcat.value = [kcat.value;kcattmp];
        kcat.ref = [kcat.ref;reftmp];
    else
        idxtmp = ismember(kcat.rxn,rxntmp);
        if kcat.value(idxtmp) == kcattmp
            kcat.ref(idxtmp) = strcat(kcat.ref(idxtmp),{'; '},reftmp);
        elseif kcat.value(idxtmp) < kcattmp
            kcat.value(idxtmp) = kcattmp;
            kcat.ref(idxtmp) = reftmp;
        end
    end
end

cd ../;
save('kcat.mat','kcat');
clear;
