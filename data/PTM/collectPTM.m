
% download from (https://downloads.thebiogrid.org/File/BioGRID/Latest-Release/BIOGRID-PTMS-LATEST.ptm.zip)
fid  = fopen('BIOGRID-PTM-4.3.194.ptmtab.txt');
data = textscan(fid,[repmat('%s ',1,17) '%s'],'Delimiter','\t','headerLines',1);
fclose(fid);
proteinlist = data{4};
PTMlist = data{10};
orglist = data{15};
idx = ismember(orglist,'Saccharomyces cerevisiae (S288c)');
proteinlist = proteinlist(idx);
PTMlist = PTMlist(idx);

list = unique(strcat(proteinlist,repmat('|',length(proteinlist),1),PTMlist));
list_split = split(list,'|');

PTMinfo = struct();
PTMinfo.protein = list_split(:,1);
PTMinfo.type = list_split(:,2);

save('PTMinfo.mat','PTMinfo');

