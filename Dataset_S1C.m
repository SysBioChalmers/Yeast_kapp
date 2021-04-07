load('kapp.mat');

maincolor = [240,59,32]/255;
heatmaplow = [255,237,160]/255;

rho2data = zeros(length(kapp.condition),length(kapp.condition));
pdata = zeros(length(kapp.condition),length(kapp.condition));
for i = 1:length(kapp.condition)
    kapp_i = kapp.values(:,i);
    idx_i = kapp_i~=0;
    rxn_i = kapp.rxn(idx_i);
    kapp_i = kapp_i(idx_i);
    for j = 1:length(kapp.condition)
        kapp_j = kapp.values(:,j);
        idx_j = kapp_j~=0;
        rxn_j = kapp.rxn(idx_j);
        kapp_j = kapp_j(idx_j);
        rxn_ij = intersect(rxn_i,rxn_j);
        [~,p] = ismember(rxn_ij,rxn_i);
        idata = kapp_i(p);
        [~,q] = ismember(rxn_ij,rxn_j);
        jdata = kapp_j(q);
        [RHOtmp,PVALtmp] = corr(log10(idata),log10(jdata),'Type','Pearson');
        rho2data(i,j) = RHOtmp^2;
        pdata(i,j) = PVALtmp;
    end
end

max(max(pdata));

figure();
xlbl = num2cell(1:1:length(kapp.condition));
xlbl_tmp = cellfun(@(x) [num2str(x),'.'],xlbl,'UniformOutput',false);
ylbl = strcat(xlbl_tmp,kapp.condition);
minclr = heatmaplow;
maxclr = maincolor;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(xlbl,ylbl,rho2data,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','k');
title(['R^2 between kapps of ' num2str(length(kapp.condition)) ' conditions']);
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[820 320 360 360]);
set(gca,'position',[0.25 0.2 0.5 0.5]);

