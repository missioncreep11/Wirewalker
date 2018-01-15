% process aquadopp with new protocole

addpath(fullfile(cd,'Toolbox/'))
addpath(fullfile(cd,'Toolbox/rsktools'))
addpath(fullfile(cd,'Toolbox/seawater'))


%% USER PART (define by user)

WWmeta.rbrpath='/Users/drew/Google Drive/NAVO/NAVO_SIO_training/rbr/';
WWmeta.aqdpath='/Users/drew/Google Drive/NAVO/NAVO_SIO_training/aqd/';
WWmeta.name_rbr='NAVO_rbr';
WWmeta.name_aqd='NAVO_aqd';
WWmeta.root_script='./';
WWmeta.WWpath='/Users/drew/Google Drive/NAVO/NAVO_SIO_training/L1/';
WWmeta.WW_name='NAVO1';
WWmeta

%% process rbr
process_rbr(WWmeta)
create_profiles_rbr(WWmeta)
create_grid_rbr(WWmeta)

%% process aqd
process_aqd_2G(WWmeta)
create_profiles_aqd_2G(WWmeta)
create_grid_aqd_2G(WWmeta)


