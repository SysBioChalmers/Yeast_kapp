load('kapp.mat');
load('kcat.mat');

maincolor = [240,59,32]/255;

row = cell(0,1);
data = zeros(0,1);
figure();
for i = 1:length(kapp.condition)
    row = [row;{[num2str(i),'.',kapp.condition{i}]}];
    kapptmp = kapp.values(:,i);
    idxtmp = kapptmp~=0;
    kapprxn = kapp.rxn(idxtmp);
    kapplist = kapptmp(idxtmp);
    rxnstmp = intersect(kcat.rxn,kapprxn);
    [~,p] = ismember(rxnstmp,kcat.rxn);
    x_kcat = kcat.value(p);
    [~,q] = ismember(rxnstmp,kapprxn);
    y_kapp = kapplist(q);
    idx_heterexp = kcat.HeterExp(p) ~= 1;
    x_kcat = x_kcat(idx_heterexp);
    y_kapp = y_kapp(idx_heterexp);
    [RHOtmp,PVALtmp] = corr(log10(x_kcat),log10(y_kapp),'Type','Pearson');
    data = [data;RHOtmp^2];
    subplot(4,7,i);
    line([-4 6],[-4 6],'Color','k');
    hold on;
    box on;
    scatter(log10(x_kcat),log10(y_kapp),3,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',maincolor,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
    xlim([-4 6]);
    ylim([-4 6]);
    xticks([]);
    yticks([]);
	text(-3.5,4.7,num2str(i),'FontSize',6,'FontName','Helvetica','FontWeight','bold');
    text(-2,-2.5,['R^2' '= ' num2str(round(RHOtmp^2,2))],'FontSize',6,'FontName','Helvetica');
end
set(gcf,'position',[400 500 360 190]);


