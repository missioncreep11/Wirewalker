function compile_RBR_telemetry(WWmeta)

deploymentpath=[WWmeta.root_data WWmeta.Cruise_name '/WW/' WWmeta.WW_name '/'];
D=load([deploymentpath 'd' num2str(1) '/L1/JohnWesleyPowell_grid.mat'],'RBRgrid');
zaxis=min([D.RBRgrid.z]):.25:max([D.RBRgrid.z]);
disp(['length zaxis' num2str(length(zaxis))])
Z=length(zaxis);
critp= 70;
critm= 20;

listfile=dir([WWmeta.telemetrypath '*.mat']);

sortdate=zeros(length(listfile),1);
for n=1:length(listfile)
    sortdate(n)=datenum(listfile(n).name(11:end-4),'mm_dd_yyyy');
end
[~,IA]=sort(sortdate);


for l=1:length(listfile)
    disp([num2str(l) 'on ' num2str(length(listfile))])
    n=IA(l);
    load([listfile(n).folder '/' listfile(n).name])
    WW_all   = structfun(@(x) x.',WW_all,'un',0);
    minP     = min(WW_all.P(:));
    WW_all.P = WW_all.P-minP;
    
    [up,~,dataup] = get_upcast_telemetry(WW_all);
    
    dup=diff(up);
    ind_prof=find(dup>1);
    RBRprofiles=struct([]);
    fields=fieldnames(dataup);
    tdata=dataup.time;
    if length(ind_prof)>2
        for i=1:length(ind_prof)-1
            for f=1:length(fields)
                wh_field=fields{f};
                if (length(tdata)==length(dataup.(wh_field)))
                    RBRprofiles{i}.(wh_field)=dataup.(wh_field)(ind_prof(i)+1:ind_prof(i+1));
                    RBRprofiles{i}.info=dataup.info;
                end
            end
        end
        %% compute N2
        disp('find depths and displacements of a few selected isopycnals')
        for i=1:length(RBRprofiles)
            if length(RBRprofiles{i}.T)>3
                RBRprofiles{i}.n2 = sw_bfrq(RBRprofiles{i}.S,...
                    RBRprofiles{i}.T,...
                    RBRprofiles{i}.P,[]);
                RBRprofiles{i}.n2=[RBRprofiles{i}.n2;RBRprofiles{i}.n2(end,:)];
            else
                RBRprofiles{i}.n2=RBRprofiles{i}.T*nan;
            end
        end
        
        Prbr=cellfun(@(x) mean(x.P),RBRprofiles);
        timerbr=cellfun(@(x) mean(x.time),RBRprofiles);
        
        timerbrOK=timerbr(Prbr>critm & Prbr<critp);
        indOK=(Prbr>critm & Prbr<critp);
        if sum(indOK)>2
            RBRprofilesOK=RBRprofiles(indOK);
            fields=fieldnames(RBRprofilesOK{1});
            for f=1:length(fields)
                wh_field=fields{f};
                if ~strcmp(wh_field,'info')
                    RBRgrid(l).(wh_field)=zeros([Z,sum(indOK)]);
                    for t=1:length(timerbrOK)
                        F=RBRprofilesOK{t}.(wh_field);
                        [Psort,I]=sort(RBRprofilesOK{t}.P,'descend');
                        P_temp=(interp_sort(Psort));
                        RBRgrid(l).(wh_field)(:,t)=interp1(P_temp,F(I),zaxis);
                    end
                end
            end
            
            %% clean n2 field
            n2=RBRgrid(l).n2;
            n2(n2<=0)=nan;
            [TT,ZZ]=meshgrid(timerbrOK,zaxis);
            TT=TT(:);ZZ=ZZ(:);
            F=scatteredInterpolant(TT(~isnan(n2)),ZZ(~isnan(n2)),n2(~isnan(n2)));
            nn2=F(TT(:),ZZ(:));
            nn2=reshape(nn2,[Z,length(timerbrOK)]);
            nn2(nn2<=0)=1e-10;
            RBRgrid(l).n2=nn2;
            
            RBRgrid(l).rho=sw_dens(RBRgrid(l).S,RBRgrid(l).T,RBRgrid(l).P);
            
            
            RBRgrid(l).time=timerbrOK;
            RBRgrid(l).z=zaxis;
            disp('coucou')
        end
    end
end

%RBRgridtot=struct([]);
fields=fieldnames(RBRgrid);
for f=1:length(fields)-3
    wh_field=fields{f};
    if ~strcmp(wh_field,'z')
        RBRgridtot.(wh_field)=[RBRgrid(:).(wh_field)];
    end
end
RBRgridtot.z=zaxis;

% %% create mask for aqd and rbr
% %  start with aqd
% [irregul_time,IA,~]=unique(RBRgridtot.time);
% dt=nanmedian(diff(irregul_time));
% timerbr=min(RBRgridtot.time):dt:max(RBRgridtot.time);
% for f=[1 2 4 5 6 8 9 10 13]
%     wh_field=fields{f};
%     RBRgridtot.(wh_field)=interp1(irregul_time,RBRgridtot.(wh_field)(:,IA).',timerbr).';
% end
% RBRgridtot.irregul_time=irregul_time;
% RBRgridtot.time=timerbr;
% 
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
datetick;
axis ij;
caxis([10 20])


save([WWmeta.compilepath '/compile_telemetry_' WWmeta.WW_name '.mat'],'RBRgridtot')

fig=gcf;
saveas(fig,'compile_RBR_telemetry.fig')



