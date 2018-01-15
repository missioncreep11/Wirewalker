function process_rbr(WWmeta)


% adding path 
addpath([WWmeta.root_script 'Toolbox']);
addpath([WWmeta.root_script 'Toolbox/rsktools']);
addpath([WWmeta.root_script 'Toolbox/gsw_matlab_v3_02']);
addpath([WWmeta.root_script 'Toolbox/gsw_matlab_v3_02/library']);
addpath([WWmeta.root_script 'Toolbox/position_scripts']);


% used in process_WW should be the name after the folder WW in dirName ;

WWmeta.rbrfile=dir(fullfile(WWmeta.rbrpath,'*.rsk'));
if length(WWmeta.rbrfile)>2;
    fprintf('Watch out \nThere is more than one rsk file\n')
    for j=1:length(filedir); disp(WWmeta.rbrfile(j).name);end
end
fprintf('read rbr file is %s\n',WWmeta.rbrfile(1).name)

disp('RSK_wrapper--- It is may take a while --- coffee time buddy?')

RSKfile= fullfile(WWmeta.rbrpath,WWmeta.rbrfile(1).name);
RSKdb=RSKopen(RSKfile);
RSKread=RSKreaddata(RSKdb);
rsk_struct_raw=RSK_struct(RSKread);
eval([WWmeta.name_rbr '=rsk_struct_raw;'])

save(fullfile(WWmeta.rbrpath,[WWmeta.name_rbr '.mat']),WWmeta.name_rbr);



