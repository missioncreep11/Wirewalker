function compile_deployement(WWmeta,range_deployment)

deploymentpath=[WWmeta.root_data WWmeta.Cruise_name '/WW/' WWmeta.WW_name '/'];
D=load([deploymentpath 'd' num2str(range_deployment(1)) '/L1/JohnWesleyPowell_grid.mat'],'RBRgrid');
fields=fieldnames(D.RBRgrid);
zaxis=min([D.RBRgrid.z]):.25:max([D.RBRgrid.z]);

for n=1:length(range_deployment)
    D=load([deploymentpath 'd' num2str(range_deployment(n)) '/L1/JohnWesleyPowell_grid.mat'],'RBRgrid');
    RBRgrid(n)=D.RBRgrid;
end

RBRgridtot=RBRgrid(1);
for f=1:length(fields)-3
    wh_field=fields{f};
    RBRgridtot.(wh_field)=[];
    for n=1:length(range_deployment)
        if strcmp(wh_field,'time')
            RBRgridtot.(wh_field)=[RBRgridtot.(wh_field) RBRgrid(n).(wh_field)];
        else
            tempo=interp1(RBRgrid(n).z,RBRgrid(n).(wh_field),zaxis);
            RBRgridtot.(wh_field)=[RBRgridtot.(wh_field) tempo];
        end
    end
end
RBRgridtot.z=zaxis;
RBRgridtot.info=RBRgrid(n).info;

% %% create mask for aqd and rbr
% %  start with aqd
% [irregul_time,IA,~]=unique(RBRgridtot.time);
% dt=nanmedian(diff(irregul_time));
% timerbr=min(RBRgridtot.time):dt:max(RBRgridtot.time);
% for f=[1 2 3 4 5 6 7 8 9 10 12 13 14]
%     wh_field=fields{f};
%     RBRgridtot.(wh_field)=interp1(irregul_time,RBRgridtot.(wh_field)(:,IA).',timerbr).';
% end
% RBRgridtot.irregul_time=irregul_time;
% RBRgridtot.time=timerbr;
% [Z,T]=size(RBRgridtot.T);
% 
% rbrmask=1+0*timerbr;
% for i=1:length(timerbr)
%     if min(abs(timerbr(i)-irregul_time))>2*dt
%         rbrmask(i)=nan;
%     end
% end
% RBRgridtot.mask=rbrmask;

figure
pcolor(RBRgridtot.time,RBRgridtot.z,RBRgridtot.T);
shading flat;
ylabel('Depth');
xlabel('time');
axis ij;
caxis([10 20])

save([WWmeta.compilepath '/compile_deployment_' WWmeta.WW_name '.mat'],'RBRgridtot')



