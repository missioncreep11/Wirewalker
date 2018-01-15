function create_profiles_aqd(WWmeta)

addpath([WWmeta.root_script 'Toolbox/seawater/'])

load([WWmeta.aqdpath WWmeta.name_aqd '.mat'])
eval(['[up,down,dataup] = get_upcast_aqd(' WWmeta.name_aqd ');'])

dup=diff(up);
ind_prof=find(dup>1);
AQDprofiles=struct([]);
fields=fieldnames(dataup);
tdata=dataup.Burst_MatlabTimeStamp;
for i=1:length(ind_prof)-1
    for f=1:length(fields)
        wh_field=fields{f};
        if (length(tdata)==length(dataup.(wh_field)))
           AQDprofiles{i}.(wh_field)=dataup.(wh_field)(ind_prof(i)+1:ind_prof(i+1));
        end
    end
end
%% add mean Corr and mean Amp field to AQDprofile
meanAmp=@(x) (nanmean([x.Burst_AmpBeam1, ...
                       x.Burst_AmpBeam2, ...
                       x.Burst_AmpBeam3, ...
                       x.Burst_AmpBeam4],2));
meanCorr=@(x) (nanmean([x.Burst_CorBeam1, ...
                        x.Burst_CorBeam2, ...
                        x.Burst_CorBeam3, ...
                        x.Burst_CorBeam4],2));
           

Corr=cellfun(meanCorr,AQDprofiles,'un',0);
Amp =cellfun(meanAmp ,AQDprofiles,'un',0);
for i=1:length(AQDprofiles)
    AQDprofiles{i}.Burst_Corr=Corr{i};
    AQDprofiles{i}.Burst_Amp =Amp{i};
end
save([WWmeta.aqdpath 'Profiles_' WWmeta.name_aqd],'AQDprofiles')

