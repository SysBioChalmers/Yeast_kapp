load('kapp.mat');
load('kcat.mat');
load('kmax.mat');
load('GEM-yeast-split.mat');
model = model_split;
clear model_split;

%% color
% maincolor = [124,81,161]/255;
% heatmaplow = [218,218,235]/255;

maincolor = [240,59,32]/255;
heatmaplow = [255,237,160]/255;

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
list_comps = cell(length(list_rxns),1);

for i = 1:length(list_rxns)
    metids = model.mets(model.S(:,ismember(model.rxns,list_rxns(i))) ~= 0);
    metcomps = cellfun(@(x) x(strfind(x,'[')+1:strfind(x,']')-1),metids,'UniformOutput',false);
    metcomps = unique(metcomps);
    if length(metcomps) == 1
        list_comps(i) = metcomps;
    else
        list_comps(i) = {'multiple comps'};
    end
end

compslist = unique(list_comps);
for i = 1:length(compslist)
    idx = ismember(list_comps,compslist(i));
    x = list_x(idx);
    y = list_y(idx);
    [RHO,PVAL] = corr(log10(x),log10(y),'Type','Pearson');
    subplot(1,length(compslist),i);
    line([-4 6],[-4 6],'Color','k');
    hold on;
    box on;
    scatter(log10(x),log10(y),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
    xlim([-4 6]);
    ylim([-4 6]);
    xticks([]);
    yticks([]);
    title(compslist{i},'FontSize',6,'FontName','Helvetica');
    text(2,0.2,['R^2' ' = ' num2str(round(RHO^2,3))],'Color','black','FontSize',6,'FontName','Helvetica');
    text(2,-0.5,['N = ' num2str(length(x))],'Color','black','FontSize',6,'FontName','Helvetica');
    text(2,-1.2,['p' ' = ' num2str(round(PVAL,7))],'Color','black','FontSize',6,'FontName','Helvetica');
end






