
load('kcat.mat');
load('kmax.mat');

rxns = intersect(kcat.rxn,kapp4.rxn);
[~,p] = ismember(rxns,kcat.rxn);
x_kcat = kcat.value(p);
[~,q] = ismember(rxns,kapp4.rxn);
y_kmax = kapp4.max(q);
idx_heterexp = kcat.HeterExp(p) == 1;
x_kcat = x_kcat(idx_heterexp);
y_kmax = y_kmax(idx_heterexp);
tablerxns = rxns(idx_heterexp);

load('PTMinfo.mat'); % check PTMs
tablePTMs = unique(PTMinfo.type)';
tableproteins = cell(length(tablerxns),1);
tableTF = false(length(tablerxns),length(tablePTMs));


for i = 1:length(tablerxns)
    prottmp = kapp4.protein(ismember(kapp4.rxn,tablerxns(i)));
    tableproteins(i) = prottmp;
    prottmp = strrep(prottmp,'( ','');
    prottmp = strrep(prottmp,' )','');
    prottmp = split(prottmp,' or ');
    for j = 1:length(tablePTMs)
        protlist = PTMinfo.protein(ismember(PTMinfo.type,tablePTMs(j)));
        tableTF(i,j) = any(ismember(prottmp,protlist));
    end
end


