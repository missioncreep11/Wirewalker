function create_grid_rbr(WWmeta)

load([WWmeta.rbrpath 'Profiles_' WWmeta.name_rbr],'RBRprofiles')


%get the normal upcast (mean P of the upcast ~ median P of all the mean P)


Prbr=cellfun(@(x) mean(x.P),RBRprofiles);
critp= max(cellfun(@(x) max(x.P),RBRprofiles))-.5*std(Prbr);
critm= min(cellfun(@(x) min(x.P),RBRprofiles))+.5*std(Prbr);


timerbr=cellfun(@(x) mean(x.time),RBRprofiles);
timerbrOK=timerbr(Prbr>critm & Prbr<critp);
indOK=(Prbr>critm & Prbr<critp);
RBRprofiles=RBRprofiles(indOK);
zaxis=0:.25:max(cellfun(@(x) max(x.P),RBRprofiles));
Z=length(zaxis);

fields=fieldnames(RBRprofiles{1});
for f=1:length(fields)
    wh_field=fields{f};
    if ~strcmp(wh_field,'info')
        RBRgrid.(wh_field)=zeros([Z,sum(indOK)]);
        for t=1:length(timerbrOK)
            F=RBRprofiles{t}.(wh_field);
            [Psort,I]=sort(RBRprofiles{t}.P,'descend');
            P_temp=(interp_sort(Psort));
            RBRgrid.(wh_field)(:,t)=interp1(P_temp,F(I),zaxis);
        end
    end
end
RBRgrid.z=zaxis;
RBRgrid.time=timerbrOK;

%% clean n2 field
n2=RBRgrid.n2;
n2(n2<=0)=nan;
[TT,ZZ]=meshgrid(RBRgrid.time,RBRgrid.z);
%TT=TT(:);ZZ=ZZ(:);
F=scatteredInterpolant(TT(~isnan(n2)),ZZ(~isnan(n2)),n2(~isnan(n2)));
nn2=F(TT(:),ZZ(:));
nn2=reshape(nn2,[length(RBRgrid.z),length(RBRgrid.time)]);
nn2(nn2<=0)=1e-10;
RBRgrid.n2=nn2;

RBRgrid.z=zaxis;
RBRgrid.time=timerbrOK;
RBRgrid.info=RBRprofiles{1}.info;

if exist([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'file')
    load([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'AQDgrid')
    save([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'RBRgrid','AQDgrid')
else
    save([WWmeta.WWpath WWmeta.WW_name '_grid.mat'],'RBRgrid')
end

