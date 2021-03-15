load('kapp.mat');
load('kcat.mat');
load('kmax.mat');

%% GO terms
% fid  = fopen('GOterm.tsv');
% data = textscan(fid,[repmat('%s ',1,4) '%s'],'Delimiter','\t','headerLines',1);
% fclose(fid);
% proteinlist = data{1};
% GOID = data{3};
% GOlist = data{4};
% clear ans data fid;

fid  = fopen('go_slim_mapping.tab');
data = textscan(fid,[repmat('%s ',1,6) '%s'],'Delimiter','\t','headerLines',1);
fclose(fid);
proteinlist = data{1};
GOlist = data{5};
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

% idx_heterexp = kcat.HeterExp(p) ~= 1;
% list_x = list_x(idx_heterexp);
% list_y = list_y(idx_heterexp);
% list_rxns = list_rxns(idx_heterexp);
% list_gpr = list_gpr(idx_heterexp);

GO_rxns = cell(0,1);
GO_x = zeros(0,1);
GO_y = zeros(0,1);
GO_GOname = cell(0,1);

for i = 1:length(list_gpr)
    gpr = list_gpr(i);
    rxn = list_rxns(i);
    
    gpr = strrep(gpr,'( ','');
    gpr = strrep(gpr,' )','');
    gpr = split(gpr,' or ');
    GOtmp = cell(0,1);
    for j = 1:length(gpr)
        if ismember(gpr(j),proteinlist)
            GOtmp = [GOtmp;unique(GOlist(ismember(proteinlist,gpr(j))))];
        end
    end
    GOtmp = unique(GOtmp);
    GO_rxns = [GO_rxns;repmat(rxn,length(GOtmp),1)];
    GO_x = [GO_x;list_x(i)*ones(length(GOtmp),1)];
    GO_y = [GO_y;list_y(i)*ones(length(GOtmp),1)];
    GO_GOname = [GO_GOname;GOtmp];
end

Group = unique(GO_GOname);
GOr = zeros(0,1);
GOp = zeros(0,1);
GOterm = cell(0,1);
GOrxns = cell(0,1);
for i = 1:length(Group)
     if sum(ismember(GO_GOname,Group(i))) >= 3
         [RHOtmp,PVALtmp] = corr(log10(GO_x(ismember(GO_GOname,Group(i)))),log10(GO_y(ismember(GO_GOname,Group(i)))),'Type','Pearson');
         if RHOtmp > 0 && PVALtmp < 0.05
             GOr = [GOr;RHOtmp];
             GOp = [GOp;PVALtmp];
             GOterm = [GOterm;Group(i)];
             GOrxns = [GOrxns;join(GO_rxns(ismember(GO_GOname,Group(i))))];
         end
     end
end





idxtest = ismember(list_rxns,split(GOrxns(5)));
scatter(log10(list_x(idxtest)),log10(list_y(idxtest)));
list_rxns(idxtest)

% for i = 1:length(list_gpr)
%     gprtmp = list_gpr{i};
%     if contains(gprtmp,' or ')
%         gprtmp = strrep(gprtmp,'( ','');
%         gprtmp = strrep(gprtmp,' )','');
%         gprtmp = strsplit(gprtmp,' or ');
%         list_gpr(i) = gprtmp(1);
%     end
% end
% list_dev = abs(log10(list_x./list_y));
% clear p q i gprtmp;
% 
% 
% idxtest = ismember(list_gpr,unique(proteinlist(ismember(GOlist,'cellular amino acid metabolic process'))));
% [RHOtmp,PVALtmp] = corr(log10(list_x(idxtest)),log10(list_y(idxtest)),'Type','Pearson');
% scatter(log10(list_x(idxtest)),log10(list_y(idxtest)));
% 
% sum(idxtest)
% 

% cutoff = 0.7; % 70%
% idx_good = (list_dev <= quantile(list_dev,cutoff));
% idx_bad = (list_dev > quantile(list_dev,cutoff));
% gpr_good = list_gpr(idx_good);
% gpr_bad = list_gpr(idx_bad);
% 
% [RHOtmp,PVALtmp] = corr(log10(list_x(idx_good)),log10(list_y(idx_good)),'Type','Pearson');
% scatter(log10(list_x(idx_good)),log10(list_y(idx_good)));


