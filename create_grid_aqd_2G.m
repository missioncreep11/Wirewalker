function create_grid_aqd_2G(WWmeta)

load([WWmeta.aqdpath 'Profiles_' WWmeta.name_aqd],'AQDprofiles')
%get the normal upcast (mean P of the upcast ~ median P of all the mean P)

Paqd=cellfun(@(x) nanmean(x.Burst_Pressure),AQDprofiles);
timeaqd=cellfun(@(x) mean(x.Burst_Time),AQDprofiles);

critp= max(cellfun(@(x) max(x.Burst_Pressure),AQDprofiles))-.5*std(Paqd);
critm= min(cellfun(@(x) min(x.Burst_Pressure),AQDprofiles))+.5*std(Paqd);
indOK=(Paqd>critm & Paqd<critp);

PaqdOK=Paqd(indOK);
timeaqdOK=timeaqd(indOK);
AQDprofilesOK=AQDprofiles(indOK);

zaxis=0:.25:max(cellfun(@(x) max(x.Burst_Pressure),AQDprofilesOK));
Z=length(zaxis);

fields=fieldnames(AQDprofiles{1});
for f=2:length(fields) % loop start at 2 because the first field is the time
    wh_field=fields{f};
    Size_field=size(double(AQDprofilesOK{1}.(wh_field)));
    AQDgrid.(wh_field)=squeeze(zeros([Z,Size_field(2),sum(indOK)]));
    for t=1:length(timeaqdOK)
        F=double(AQDprofilesOK{t}.(wh_field));
        if length(size(F))==3;F=F(:,:,1);end
        [Psort,~]=sort(AQDprofilesOK{t}.Burst_Pressure,'descend');
        [P_temp,IA,~]=unique(Psort);
        if size(F,2)==1
            AQDgrid.(wh_field)(:,t)=interp1(P_temp,medfilt1(F(IA,:),10),zaxis);
        else
            AQDgrid.(wh_field)(:,:,t)=interp1(P_temp,medfilt1(F(IA,:),10),zaxis);
        end
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


