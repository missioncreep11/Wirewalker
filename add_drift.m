function add_drift(WWmeta)

%% add positions of WW_name
    fprintf('Add WW trajectory\n')
    GPS=load(WWmeta.GPS);
    load([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'AQDgrid')
    
    WW_pos(:,1)=GPS.([WWmeta.WW_name '_pos']).GPSLongitude;
    WW_pos(:,2)=GPS.([WWmeta.WW_name '_pos']).GPSLatitude;
    WW_pos(:,3)=GPS.([WWmeta.WW_name '_pos']).DataDateTime;

    lon=WW_pos(GPS.([WWmeta.WW_name '_pos']).gpsQuality==3,1);
    lat=WW_pos(GPS.([WWmeta.WW_name '_pos']).gpsQuality==3,2);
    time1=WW_pos(GPS.([WWmeta.WW_name '_pos']).gpsQuality==3,3);

    AQDgrid.lat=interp1(time1,lat,AQDgrid.time)';
    AQDgrid.lon=interp1(time1,lon,AQDgrid.time)';


    [dist1,az]=m_lldist_az(AQDgrid.lon, AQDgrid.lat);
    dist1=[0 dist1'];
    az=[0 az'];
    dist1=cumsum(dist1);

    AQDgrid.dist=dist1;
    AQDgrid.drift=nan(size(AQDgrid.dist));
    AQDgrid.drift(2:end)=diff(AQDgrid.dist*1000)./(diff(AQDgrid.time)*86400);
    AQDgrid.drift(1)=AQDgrid.drift(2);
    AQDgrid.n_drift=AQDgrid.drift.*cosd(az);
    AQDgrid.e_drift=AQDgrid.drift.*sind(az);

    AQDgrid.n_drift=conv2(AQDgrid.n_drift,gausswin(12)'./sum(gausswin(12)),'same');
    AQDgrid.e_drift=conv2(AQDgrid.e_drift,gausswin(12)'./sum(gausswin(12)),'same');
    
    AQDgrid.e_abs=AQDgrid.Burst_VelEast+repmat(AQDgrid.e_drift,length(AQDgrid.z),1);
    AQDgrid.n_abs=AQDgrid.Burst_VelNorth+repmat(AQDgrid.n_drift,length(AQDgrid.z),1);
    
    if exist([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'file')
        load([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'RBRgrid')
        save([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'RBRgrid','AQDgrid')
    else
        save([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'AQDgrid')
    end
end
    
