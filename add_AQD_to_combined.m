function add_AQD_to_combined(WWmeta,range_deployment)


deploymentpath=[WWmeta.root_data WWmeta.Cruise_name '/WW/' WWmeta.WW_name '/'];

load(fullfile(deploymentpath,'compile_deployment',['compile_',WWmeta.WW_name,'.mat']))
if ~isfield(RBRgrid,'u')
    RBRgrid.u   = NaN.*RBRgrid.T;
    RBRgrid.v   = NaN.*RBRgrid.T;
    RBRgrid.w   = NaN.*RBRgrid.T;
    RBRgrid.amp = NaN.*RBRgrid.T;
end

for ii =range_deployment
    load(fullfile(deploymentpath,['d' num2str(ii)], 'L1/JohnWesleyPowell_grid.mat'),'AQDgrid');
    try
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
            id = findnearest(time(jj),RBRgrid.time);
            RBRgrid.u(:,id)     = interp1(z,u(:,jj),RBRgrid.z);
            RBRgrid.v(:,id)     = interp1(z,v(:,jj),RBRgrid.z);
            RBRgrid.w(:,id)     = interp1(z,w(:,jj),RBRgrid.z);
            RBRgrid.amp(:,id)   = interp1(z,amp(:,jj),RBRgrid.z);
        end
        
    catch ME
        disp(['AQD not incorporated for deployment ',num2str(ii)])
    end
end

save(fullfile(deploymentpath,'compile_deployment',['compile_',WWmeta.WW_name,'.mat']),'RBRgrid')
