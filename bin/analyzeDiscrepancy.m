load('kapp.mat');
load('kcat.mat');
load('kmax.mat');
load('GEM-yeast.mat');

%% GO terms
fid  = fopen('GOterm.tsv');
data = textscan(fid,[repmat('%s ',1,4) '%s'],'Delimiter','\t','headerLines',1);
fclose(fid);
proteinlist = data{1};
GOlist = data{4};
pathwaylist = data{5};
clear ans data fid;

%% color
maincolor = [124,81,161]/255;
heatmaplow = [218,218,235]/255;

%% kmax vs kcat

% kapp4.rxn = kapp4.rxn(~contains(kapp4.protein,'or'));
% kapp4.max = kapp4.max(~contains(kapp4.protein,'or'));

kcat_rxn = kcat.rxn;
kcat_num = kcat.value;

kapp_num = kapp4.max;
kapp_gpr = kapp4.protein;
kapp_rxn = kapp4.rxn;


% i = 1;
% kapptmp = kapp.values(:,i);
% idxtmp = kapptmp~=0;
% kapp_rxn = kapp.rxn(idxtmp);
% kapp_num = kapptmp(idxtmp);
% kapp_gpr = kapp.protein(idxtmp);



list_rxns = intersect(kcat_rxn,kapp_rxn);
[~,p] = ismember(list_rxns,kcat_rxn);
x = kcat_num(p);
[~,q] = ismember(list_rxns,kapp_rxn);
y = kapp_num(q);
list_gpr = kapp_gpr(q);
list_dev = abs(log10(x./y));
clear p q;

GO.rxns = cell(0,1);
GO.dev = zeros(0,1);
GO.GO = cell(0,1);

PW.rxns = cell(0,1);
PW.dev = zeros(0,1);
PW.pathway = cell(0,1);

for i = 1:length(list_gpr)
    gpr = list_gpr(i);
    rxn = list_rxns(i);
    dev = list_dev(i);
    
    gpr = strrep(gpr,'( ','');
    gpr = strrep(gpr,' )','');
    gpr = split(gpr,' or ');
    GOtmp = cell(0,1);
    pathwaytmp = cell(0,1);
    for j = 1:length(gpr)
        if ismember(gpr(j),proteinlist)
            GOtmp = [GOtmp;unique(GOlist(ismember(proteinlist,gpr(j))))];
            pathwaytmp = [pathwaytmp;unique(pathwaylist(ismember(proteinlist,gpr(j))))];
        end
    end
    GOtmp = unique(GOtmp);
    pathwaytmp = unique(pathwaytmp);
    GO.rxns = [GO.rxns;repmat(rxn,length(GOtmp),1)];
    GO.dev = [GO.dev;dev*ones(length(GOtmp),1)];
    GO.GO = [GO.GO;GOtmp];
    PW.rxns = [PW.rxns;repmat(rxn,length(pathwaytmp),1)];
    PW.dev = [PW.dev;dev*ones(length(pathwaytmp),1)];
    PW.pathway = [PW.pathway;pathwaytmp];
    
end

Group = unique(GO.GO);
GOdev = zeros(0,1);
GOterm = cell(0,1);
for i = 1:length(Group)
     if sum(ismember(GO.GO,Group(i))) >= 10
         GOdev = [GOdev;mean(GO.dev(ismember(GO.GO,Group(i))))];
         GOterm = [GOterm;Group(i)];
     end
end



rxntest = GO.rxns(ismember(GO.GO,GOterm(GOdev == max(GOdev))));
[~,b] = ismember(rxntest,list_rxns);
[RHO,PVAL] = corr(log10(x(b)),log10(y(b)),'Type','Pearson');
scatter(log10(x(b)),log10(y(b)));


Group = unique(PW.pathway);
pathwaydev = zeros(0,1);
pathway = cell(0,1);
for i = 1:length(Group)
     if sum(ismember(PW.pathway,Group(i))) >= 5
         pathwaydev = [pathwaydev;mean(PW.dev(ismember(PW.pathway,Group(i))))];
         pathway = [pathway;Group(i)];
     end
end








