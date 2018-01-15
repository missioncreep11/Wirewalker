function WWgrid = combine_addfields_WW_deployments(WWmeta,deployment_range);
% function WWgrid = combine_addfields_WW_deployments(WWmeta,deployment_range);
% This function combines data for each deployment, including RBR data from
% .rsk files, RBR data from the telemetered data files where needed, and
% Aquadopp data if it has been processed. 
% output WWgrid files contain the fields:
%     'T'
%     'C'
%     'DO'
%     'BScat'
%     'F_chla'
%     'F_CDOM'
%     'time'
%     'z'
%     'u'
%     'v'
%     'w'
%     'amp'
%     'L_ot'
%     'Lmin'
%     'eps_ot'
%     'n2'
%     'rho'
%     'S'
% These files are saved for each deployement at the top level of the
% deployment directory (.../d##/WWgrid.mat)
    
load([WWmeta.compilepath '/compile_telemetry_' WWmeta.WW_name '.mat'],'RBRgridtot');
[RBRgridtot.time,id] = sort(RBRgridtot.time); id = id(find(unique(RBRgridtot.time)));
RBRgridtot.DO = RBRgridtot.DO(:,id);
RBRgridtot.BScat = RBRgridtot.BScat(:,id);
RBRgridtot.F_chla = RBRgridtot.F_chla(:,id);
RBRgridtot.F_CDOM = RBRgridtot.F_CDOM(:,id);
RBRgridtot.C = RBRgridtot.C(:,id);
RBRgridtot.T = RBRgridtot.T(:,id);
lz = length(RBRgridtot.z);
lt = length(RBRgridtot.time);

load(fullfile(WWmeta.data_path,'Index.mat'))
if nargin<2
    deployment_range = 1:length(Index.start);
end

%%
for ii = deployment_range

    %%
    depname = ['d',num2str(ii)];
    load(fullfile(WWmeta.data_path,depname,'L1',[WWmeta.WW_name,'_grid.mat']))

    WWgrid.T = RBRgrid.T;
    WWgrid.C = RBRgrid.C;
    %     WWgrid.rho = RBRgrid.rho; %we will despike and add these fields later
    %     WWgrid.S = RBRgrid.v7;
    WWgrid.DO = RBRgrid.v1;
    WWgrid.BScat = RBRgrid.v2;
    WWgrid.F_chla = RBRgrid.v3;
    WWgrid.F_CDOM = RBRgrid.v4;
    WWgrid.time = RBRgrid.time;
    WWgrid.z = RBRgrid.z';
    lzd = length(WWgrid.z);
    ltd = length(WWgrid.time);

    clear RBRgrid

    %% Telemetry
    % see if there are profiles in the telemetry data between the end of
    % the prior deployment and the beginning of the next RBR grid
    if ii>1
        id = find(RBRgridtot.time>Index.end(ii-1) & RBRgridtot.time<WWgrid.time(1));
        lt = length(id);
        if ~isempty(id)
            if lz>lzd;
                WWgrid.DO       = [RBRgridtot.DO(:,id)      [WWgrid.DO;     NaN*ones(lz-lzd,ltd)]];
                WWgrid.BScat    = [RBRgridtot.BScat(:,id)   [WWgrid.BScat;  NaN*ones(lz-lzd,ltd)]];
                WWgrid.F_chla   = [RBRgridtot.F_chla(:,id)  [WWgrid.F_chla; NaN*ones(lz-lzd,ltd)]];
                WWgrid.F_CDOM   = [RBRgridtot.F_CDOM(:,id)  [WWgrid.F_CDOM; NaN*ones(lz-lzd,ltd)]];
                WWgrid.T        = [RBRgridtot.T(:,id)       [WWgrid.T;      NaN*ones(lz-lzd,ltd)]];
                WWgrid.C        = [RBRgridtot.C(:,id)       [WWgrid.C;      NaN*ones(lz-lzd,ltd)]];
                WWgrid.z        = RBRgridtot.z';


            elseif lz<lzd
                WWgrid.DO       = [[RBRgridtot.DO(:,id);        NaN*ones(lzd-lz,lt)]    WWgrid.DO     ];
                WWgrid.BScat    = [[RBRgridtot.BScat(:,id);     NaN*ones(lzd-lz,lt)]    WWgrid.BScat  ];
                WWgrid.F_chla   = [[RBRgridtot.F_chla(:,id);    NaN*ones(lzd-lz,lt)]    WWgrid.F_chla ];
                WWgrid.F_CDOM   = [[RBRgridtot.F_CDOM(:,id);    NaN*ones(lzd-lz,lt)]    WWgrid.F_CDOM ];
                WWgrid.T        = [[RBRgridtot.T(:,id);         NaN*ones(lzd-lz,lt)]    WWgrid.T      ];
                WWgrid.C        = [[RBRgridtot.C(:,id);         NaN*ones(lzd-lz,lt)]    WWgrid.C      ];
            else
                WWgrid.DO       = [RBRgridtot.DO(:,id)     WWgrid.DO     ];
                WWgrid.BScat    = [RBRgridtot.BScat(:,id)  WWgrid.BScat  ];
                WWgrid.F_chla   = [RBRgridtot.F_chla(:,id)	WWgrid.F_chla ];
                WWgrid.F_CDOM   = [RBRgridtot.F_CDOM(:,id)	WWgrid.F_CDOM ];
                WWgrid.T        = [RBRgridtot.T(:,id)      WWgrid.T      ];
                WWgrid.C        = [RBRgridtot.C(:,id)      WWgrid.C      ];
            end

            [WWgrid.time, id]     = unique([RBRgridtot.time(id) WWgrid.time]);
            WWgrid.DO       = WWgrid.DO(:,id);
            WWgrid.BScat    = WWgrid.BScat(:,id);
            WWgrid.F_chla   = WWgrid.F_chla(:,id);
            WWgrid.F_CDOM   = WWgrid.F_CDOM(:,id);
            WWgrid.T        = WWgrid.T(:,id);
            WWgrid.C        = WWgrid.C(:,id);

        end
    end
    clear lzd ltd id

    %% Add velocities

    WWgrid.u = NaN.*WWgrid.T;
    WWgrid.v = NaN.*WWgrid.T;
    WWgrid.w = NaN.*WWgrid.T;
    WWgrid.amp = NaN.*WWgrid.T;

    if exist('AQDgrid','var')
        id = find(AQDgrid.Burst_Corr>30);
        mask = NaN.*AQDgrid.Burst_Amp; mask(id) = 1;

        u = AQDgrid.Burst_VelEast.*mask;
        v = AQDgrid.Burst_VelNorth.*mask;
        try
            w = AQDgrid.Burst_VelUp.*mask;
        catch ME
            try
                w = AQDgrid.Burst_Up.*mask;
            catch ME
                disp('no w')
            end
        end
        amp = AQDgrid.Burst_Amp;
        time = AQDgrid.time;
        z = AQDgrid.z;

        clear AQDgrid

        for jj = 1:length(time)
            id = findnearest(time(jj),WWgrid.time);
            WWgrid.u(:,id)     = interp1(z,u(:,jj),WWgrid.z');
            WWgrid.v(:,id)     = interp1(z,v(:,jj),WWgrid.z');
            WWgrid.w(:,id)     = interp1(z,w(:,jj),WWgrid.z');
            WWgrid.amp(:,id)   = interp1(z,amp(:,jj),WWgrid.z');
        end
    else
        disp(['AquaDopp for deployment ',num2str(ii),' not processed.'])
    end
    clear u v w amp time z id mask ME

    %% add Thorpe Scale overturns and N2

    load(fullfile(WWmeta.data_path,depname,'rbr',[WWmeta.WW_name,'_rbr_',depname,'.mat']))
    eval(['RBRgrid = ',WWmeta.WW_name,'_rbr_',depname,'; clear ', WWmeta.WW_name,'_rbr_',depname]);

%%
P = fastsmooth(RBRgrid.P,30,5,1);
idup = find(P>2 & diffs(P)<0);
pid = find(diff(idup)>1000); pid = [1;pid+1;length(idup)];
T = RBRgrid.T(idup);
P = RBRgrid.P(idup);
C = RBRgrid.C(idup);
time = RBRgrid.time(idup);
%     clf; plot(RBRgrid.time,RBRgrid.P,time,P,WWgrid.time,50.*ones(size(WWgrid.time)),'o');
%     hold on; plot(RBRgrid.time(idup(pid)), RBRgrid.P(idup(pid)),'ok')

%%


WWgrid.L_ot = NaN.*WWgrid.T;
WWgrid.Lmin = NaN.*WWgrid.T;
WWgrid.eps_ot = NaN.*WWgrid.T;
WWgrid.n2 = NaN.*WWgrid.T;
WWgrid.rho = NaN.*WWgrid.T;
WWgrid.S = NaN.*WWgrid.T;

for jj = 1:length(WWgrid.time)
    
    id = findnearest(WWgrid.time(jj),time(pid),-1);
    if ~isempty(id)
        id = fliplr(pid(id):pid(id+1));
        Ttmp = T(id)';
        Ctmp = C(id)';
        Ptmp = P(id)';
    else
        if ~exist('WW_all','var') || WW_all.time(end)<WWgrid.time(jj)
            load(fullfile(WWmeta.telemetrypath,[WWmeta.sn,'_',datestr(WWgrid.time(jj),'mm'),'_',...
                num2str(day(WWgrid.time(jj))),'_',datestr(WWgrid.time(jj),'yyyy'),'.mat']))
            P2 = fastsmooth(WW_all.Ps,30,5,1)';
            idup2 = find(P2>2 & diffs(P2)<0);
            pid2 = find(diff(idup2)>1000); pid2 = [1;pid2+1;length(idup2)];
            T2 = WW_all.T(idup2);
            P2 = WW_all.Ps(idup2);
            C2 = WW_all.C(idup2);
            time2 = WW_all.time(idup2);
        end
        id2 = findnearest(WWgrid.time(jj),time2(pid2),-1);
        id2 = fliplr(pid2(id2):pid2(id2+1));
        Ttmp = T2(id2);
        Ctmp = C2(id2);
        Ptmp = P2(id2);
    end
    
    try
        [data, ~,~,~] = ww_despike(Ttmp,Ctmp,Ptmp,0);
        
        Ttmp = data.T;
        Stmp = data.S;
        Ptmp = data.P;
        
        WWgrid.S(:,jj) = interp1(Ptmp,Stmp,WWgrid.z);
        WWgrid.T(:,jj) = interp1(Ptmp,Ttmp,WWgrid.z);
        
        [Epsout,Lmin,Lot,~,~,n2,pout]=compute_overturns_WW(Ptmp',Ttmp',Stmp',...
            'lat',32,'lon',-117,'usetemp',0,'minotsize',1,'sigma',5e-6,'runlmin',sqrt(6));
        
        WWgrid.L_ot(:,jj) = interp1(Ptmp,Lot,WWgrid.z);
        WWgrid.Lmin(:,jj) = interp1(Ptmp,Lmin,WWgrid.z);
        WWgrid.eps_ot(:,jj) = interp1(Ptmp,Epsout,WWgrid.z);
        WWgrid.n2(:,jj) = interp1(pout,n2,WWgrid.z);
        WWgrid.rho(:,jj) = sw_pden(WWgrid.S(:,jj)',WWgrid.T(:,jj)',WWgrid.z',0);
        
    catch ME
        disp(['Profile failed for file ',num2str(ii),' profile ',num2str(jj)]);
    end
    clear id id2 Ttmp Stmp Ctmp Ptmp data Epsout Lot Lmin n2 pout
    
    
    
end

clear S T C P S2 T2 C2 P2 pid idup pid2 idup2
save(fullfile(WWmeta.data_path,depname,'WWgrid.mat'),'WWgrid')

Index.fullgridpath{ii}=fullfile(WWmeta.data_path,depname,'WWgrid.mat');
Index.start(ii) = WWgrid.time(1);
Index.nprofiles(ii) = length(WWgrid.time);
save(fullfile(WWmeta.data_path,'Index.mat'),'Index')

clear WWgrid



end

