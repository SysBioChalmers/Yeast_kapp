load('GEM-yeast-split.mat');
rxnlist = model_split.rxns;
grlist = model_split.grRules;
rowlist = cell(1,0);
valuelist = zeros(length(rxnlist),0);
fluxlist = zeros(length(rxnlist),0);
abundlist = zeros(length(rxnlist),0);

% import molecular weight
[num,txt,~] = xlsread('UniProt.xlsx');
geneUniProt = txt(2:end,1);
MWUniProt = num;
clear num txt;

cd Fluxes/; 
file = dir('*.mat');
for i = 1:length(file)
    display([num2str(i) '/' num2str(length(file))]);
    filename = file(i).name;
    load(filename);
    if Fluxes.deviation < 0.1
        condid = strrep(filename,'Fluxes_','');
        condid = strrep(condid,'.mat','');
        rowlist = [rowlist {condid}];
        
        sheetName = condid(1:strfind(condid,'_')-1);
        [num,txt,~] = xlsread('ProteomicsFlux.xlsx',sheetName);
        % remove proteins with extremely low abundance (< 1 percentile)
        numtmp = num(:);
        numtmp = numtmp(numtmp ~= 0);
        num(num < quantile(numtmp,0.01)) = 0;
        
        head = txt(1,2:end);
        proteinList = txt(2:end,1);
        abundList = num(:,ismember(head,condid));
        
        fluxes = Fluxes.pFBA;
        fluxes(abs(fluxes)<1e-5) = 0; %﻿the absolute flux value should surpass 0.00001
        kapplist_tmp = zeros(length(rxnlist),1);
        fluxeslist_tmp = zeros(length(rxnlist),1);
        abundlist_tmp = zeros(length(rxnlist),1);
        for j = 1:length(rxnlist)
            flux_tmp = fluxes(ismember(Fluxes.model.rxns,rxnlist(j)));
            
            if flux_tmp > 0 && ~ismember(grlist(j),'') && ~contains(grlist(j),' and ')
                tot_flux_tmp = Fluxes.pFBA(ismember(Fluxes.model.rxns,'r_tot_flux'));
                model_tmp = Fluxes.model;
                model_tmp = changeRxnBounds(model_tmp,'r_tot_flux',tot_flux_tmp*0.9999,'l');
                model_tmp = changeRxnBounds(model_tmp,'r_tot_flux',tot_flux_tmp*1.0001,'u');
                model_tmp = changeObjective(model_tmp,rxnlist{j});
                solmin = optimizeCbModel(model_tmp,'min');
                
                % remove data when min of FVA is zero
                if solmin.f > 0
                    gr_tmp = grlist{j};
                    gr_tmp = strrep(gr_tmp,'(','');
                    gr_tmp = strrep(gr_tmp,')','');
                    gr_tmp = strrep(gr_tmp,' or ',' ');
                    gr_tmp = strtrim(gr_tmp);
                    gr_tmp = strsplit(gr_tmp);
                    gr_tmp = unique(gr_tmp)';

                    [a,~] = ismember(gr_tmp,proteinList);
                    if any(a)
                        gr_tmp = gr_tmp(a);
                        if strcmp(sheetName,'Lahtvee2017')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            mmol_gCDW = sum(abund_tmp*1e12/6.02e23)*1000;
                            kapp_tmp = flux_tmp/mmol_gCDW/3600; %/s
                        elseif strcmp(sheetName,'Yu2020')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            mmol_gCDW = sum(abund_tmp*1e3/1e15)*1000;
                            kapp_tmp = flux_tmp/mmol_gCDW/3600; %/s
                        elseif strcmp(sheetName,'DiBartolomeo2020')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            [~,mw_idx] = ismember(gr_tmp,geneUniProt);
                            mw_tmp = MWUniProt(mw_idx);
                            mmol_gCDW = sum(abund_tmp./mw_tmp)*1000;
                            kapp_tmp = flux_tmp/mmol_gCDW/3600; %/s
                        elseif strcmp(sheetName,'Yu2021')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            mmol_gCDW = sum(abund_tmp*1e3/1e15)*1000;
                            kapp_tmp = flux_tmp/mmol_gCDW/3600; %/s
                        end
                        kapplist_tmp(j,1) = kapp_tmp;
                        abundlist_tmp(j,1) = mmol_gCDW;
                        fluxeslist_tmp(j,1) = flux_tmp;
                    end
                end
            end
        end
        valuelist = [valuelist kapplist_tmp];
        fluxlist = [fluxlist fluxeslist_tmp];
        abundlist = [abundlist abundlist_tmp];
    end
end
cd ../;

idx1 = any(valuelist,2); % remain rxns that have kapp at >= 1 condition
rxnlist = rxnlist(idx1);
grlist = grlist(idx1);
valuelist = valuelist(idx1,:);
valuelist(valuelist == inf) = 0;
fluxlist = fluxlist(idx1,:);
abundlist = abundlist(idx1,:);

kapp_raw = struct();
kapp_raw.values = valuelist;
kapp_raw.rxn = rxnlist;
kapp_raw.protein = grlist;
kapp_raw.condition = rowlist;
kapp_raw.fluxes = fluxlist;
kapp_raw.protein_conc = abundlist;

cd ../;
save('kapp_raw.mat','kapp_raw');
clear;




