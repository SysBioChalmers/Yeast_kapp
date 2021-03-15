
function fluxes = simulateCarbonSourceFixMu(model,CSexlist,maxmulist)

% set YNB with amino acids
list =  {'r_1893_REV'; ... % L-histidine exchange
    'r_1902_REV'; ... % L-methionine exchange
    'r_1912_REV'; ... % L-tryptophan exchange
    'r_1879_REV'; ... % L-arginine exchange
    'r_1881_REV'; ... % L-aspartate exchange
    'r_1897_REV'; ... % L-isoleucine exchange
    'r_1899_REV'; ... % L-leucine exchange
    'r_1900_REV'; ... % L-lysine exchange
    'r_1903_REV'; ... % L-phenylalanine exchange
    'r_1911_REV'; ... % L-threonine exchange
    'r_1913_REV'; ... % L-tyrosine exchange
    'r_1914_REV'};    % L-valine exchange

% GECKO used 2 as upper bounds. Here used 0.04 so that the predicted max
% growth rate is close to experimental value.
model = changeRxnBounds(model,list,ones(length(list),1)*0.04,'u');

model = changeRxnBounds(model,'r_1714_REV',0,'u');
model = changeRxnBounds(model,'prot_pool_exchange',1000,'u');
model = changeRxnBounds(model,'prot_pool_exchange',0,'l');

model = changeObjective(model,'prot_pool_exchange');

fluxes = zeros(length(model.rxns),length(CSexlist));
for i = 1:length(CSexlist)
    modeltmp = model;
    cstmp = CSexlist{i};
    modeltmp = changeRxnBounds(modeltmp,cstmp,1000,'u');
    modeltmp = changeRxnBounds(modeltmp,'r_2111',maxmulist(i),'b');
    soltmp = optimizeCbModel(modeltmp,'min');
    if ~isempty(soltmp.x)
        fluxes(:,i) = soltmp.x;
    end
end
