function create_grid_aqd(WWmeta)

load([WWmeta.aqdpath 'Profiles_' WWmeta.name_aqd],'AQDprofiles')
%get the normal upcast (mean P of the upcast ~ median P of all the mean P)

Paqd=cellfun(@(x) nanmean(x.Burst_Pressure),AQDprofiles);
timeaqd=cellfun(@(x) mean(x.Burst_MatlabTimeStamp),AQDprofiles);

critp= max(cellfun(@(x) max(x.Burst_Pressure),AQDprofiles))-.5*std(Paqd);
critm= min(cellfun(@(x) min(x.Burst_Pressure),AQDprofiles))+.5*std(Paqd);
indOK=(Paqd>critm & Paqd<critp);

PaqdOK=Paqd(indOK);
timeaqdOK=timeaqd(indOK);
AQDprofilesOK=AQDprofiles(indOK);

zaxis=0:.25:max(cellfun(@(x) max(x.Burst_Pressure),AQDprofilesOK));
Z=length(zaxis);

fields=fieldnames(AQDprofiles{1});
for f=1:length(fields)
    wh_field=fields{f};
    AQDgrid.(wh_field)=zeros([Z,sum(indOK)]);
    for t=1:length(timeaqdOK)
        F=AQDprofilesOK{t}.(wh_field);
        [Psort,I]=sort(AQDprofilesOK{t}.Burst_Pressure,'descend');
        P_temp=(interp_sort(Psort));
        AQDgrid.(wh_field)(:,t)=interp1(P_temp,medfilt1(F(I),10),zaxis);
    end
end
AQDgrid.z=zaxis;
AQDgrid.time=timeaqdOK;
if exist([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'file')
    load([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'RBRgrid')
    save([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'RBRgrid','AQDgrid')
else
    save([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'AQDgrid')
end


