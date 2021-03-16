load('kapp.mat');
load('kcat.mat');
load('kmax.mat');

%% pathway

fid  = fopen('kegg_pathway_yeast.txt');
data = textscan(fid,'%s %s','Delimiter','\t','headerLines',0);
fclose(fid);
pathwayid = data{1};
pathwayname = data{2};
pathwayname = strrep(pathwayname,' - Saccharomyces cerevisiae (budding yeast)','');
clear ans data fid;

fid  = fopen('kegg_gene_yeast.txt');
data = textscan(fid,'%s %s','Delimiter','\t','headerLines',0);
fclose(fid);
genelist = data{1};
genelist = strrep(genelist,'sce:','');
pathwayidlist = data{2};
for i = 1:length(pathwayidlist)
    pathwayidlist(i) = pathwayname(ismember(pathwayid,pathwayidlist(i)));
end
clear ans data fid i;

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

% i = 4;
% kapptmp = kapp.values(:,i);
% idxtmp = kapptmp~=0;
% kapp_rxn = kapp.rxn(idxtmp);
% kapp_num = kapptmp(idxtmp);
% kapp_gpr = kapp.protein(idxtmp);

list_rxns = intersect(kcat_rxn,kapp_rxn);
[~,p] = ismember(list_rxns,kcat_rxn);
list_x = kcat_num(p);
[~,q] = ismember(list_rxns,kapp_rxn);
list_y = kapp_num(q);
list_gpr = kapp_gpr(q);

idx_heterexp = kcat.HeterExp(p) ~= 1;
list_x = list_x(idx_heterexp);
list_y = list_y(idx_heterexp);
list_rxns = list_rxns(idx_heterexp);
list_gpr = list_gpr(idx_heterexp);

pw_rxns = cell(0,1);
pw_x = zeros(0,1);
pw_y = zeros(0,1);
pw_pathwayname = cell(0,1);

for i = 1:length(list_gpr)
    gpr = list_gpr(i);
    rxn = list_rxns(i);
    
    gpr = strrep(gpr,'( ','');
    gpr = strrep(gpr,' )','');
    gpr = split(gpr,' or ');
    pwtmp = cell(0,1);
    for j = 1:length(gpr)
        if ismember(gpr(j),genelist)
            pwtmp = [pwtmp;unique(pathwayidlist(ismember(genelist,gpr(j))))];
        end
    end
    pwtmp = unique(pwtmp);
    pw_rxns = [pw_rxns;repmat(rxn,length(pwtmp),1)];
    pw_x = [pw_x;list_x(i)*ones(length(pwtmp),1)];
    pw_y = [pw_y;list_y(i)*ones(length(pwtmp),1)];
    pw_pathwayname = [pw_pathwayname;pwtmp];
end

Group = unique(pw_pathwayname);
PWr = zeros(0,1);
PWp = zeros(0,1);
Pathway = cell(0,1);
PWrxns = cell(0,1);
for i = 1:length(Group)
     if sum(ismember(pw_pathwayname,Group(i))) >= 5
         [RHOtmp,PVALtmp] = corr(log10(pw_x(ismember(pw_pathwayname,Group(i)))),log10(pw_y(ismember(pw_pathwayname,Group(i)))),'Type','Pearson');
         if RHOtmp > 0 && PVALtmp < 0.05
             PWr = [PWr;RHOtmp];
             PWp = [PWp;PVALtmp];
             Pathway = [Pathway;Group(i)];
             PWrxns = [PWrxns;join(pw_rxns(ismember(pw_pathwayname,Group(i))))];
         end
     end
end

idxtest = ismember(list_rxns,split(PWrxns(2)));
scatter(log10(list_x(idxtest)),log10(list_y(idxtest)));
list_rxns(idxtest)



