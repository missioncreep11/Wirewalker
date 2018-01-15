% process aquadopp with new protocole

addpath(fullfile(cd,'Toolbox/'))
addpath(fullfile(cd,'Toolbox/rsktools'))
addpath(fullfile(cd,'Toolbox/seawater'))


%% USER PART (define by user)
WWmeta.root_data='/Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/';
WWmeta.root_script='/Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/WireWalker_master/';
WWmeta.Cruise_name='NAVO'; % 
WWmeta.WW_name='NAVO'; % 
WWmeta.deployement='d1';

%% create path
WWmeta.WWpath=sprintf('%s/%s/WW/%s/%s/L1/',WWmeta.root_data,...
    WWmeta.Cruise_name,...
    WWmeta.WW_name,...
    WWmeta.deployement);

%%
WWmeta.aqdpath=sprintf('%s%s/WW/%s/%s/aqd/',WWmeta.root_data,...
    WWmeta.Cruise_name,...
    WWmeta.WW_name,...
    WWmeta.deployement);
WWmeta.name_aqd=[WWmeta.WW_name '_aqd_' WWmeta.deployement];

WWmeta.rbrpath=sprintf('%s%s/WW/%s/%s/rbr/',WWmeta.root_data,...
    WWmeta.Cruise_name,...
    WWmeta.WW_name,...
    WWmeta.deployement);
WWmeta.name_rbr=[WWmeta.WW_name '_rbr_' WWmeta.deployement];

WWmeta.figure_path=[WWmeta.root_data 'FIGURES/'];


%% process rbr
process_rbr(WWmeta)
create_profiles_rbr(WWmeta)
create_grid_rbr(WWmeta)

%% process aqd
process_aqd_2G(WWmeta)
create_profiles_aqd_2G(WWmeta)
create_grid_aqd_2G(WWmeta)


