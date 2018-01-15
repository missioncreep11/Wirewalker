% merged telemetry and JohnWesleyPowell

RBR1=load([WWmeta.compilepath '/compile_telemetry_' WWmeta.WW_name '.mat'],'RBRgridtot');
RBR2=load([WWmeta.compilepath '/compile_deployment_' WWmeta.WW_name '.mat'],'RBRgridtot');

% RBR1 index for field with different name for RBR1 and RBR2
indRBR1.C      =  1; %RBR2 C
indRBR1.T      =  2; %RBR2 T
indRBR1.P      =  3; %RBR2 C
indRBR1.DO     =  4; %RBR2 v1
indRBR1.BScat  =  5; %RBR2 v2  
indRBR1.F_chla =  6; %RBR2 v3
indRBR1.F_CDOM =  7; %RBR2 v4
indRBR1.Ps     =  8; %RBR2 v5
indRBR1.S      = 13; %RBR2 S


[time,IA,IB]=unique([RBR1.RBRgridtot.time RBR2.RBRgridtot.time]);


fields=fieldnames(indRBR1);
fields2=fieldnames(RBR2.RBRgridtot);
RBRgrid.time=time;
RBRgrid.z=RBR1.RBRgridtot.z;

for f=1:length(fields)
    wh_field=fields{f};
    wh_field2=fields2{indRBR1.(wh_field)};
    tempo=[RBR1.RBRgridtot.(wh_field) RBR2.RBRgridtot.(wh_field2)];
    RBRgrid.(wh_field)=tempo(:,IA);
end

save([WWmeta.compilepath '/compile_' WWmeta.WW_name '.mat'],'RBRgrid')


[Z,T]=size(RBRgrid.T);

ax(1)=subplot(311);
pcolor(RBR1.RBRgridtot.time,RBR1.RBRgridtot.z,RBR1.RBRgridtot.T);
shading flat 
axis ij
ylabel('depth')
xlabel('time')
datetick

ax(2)=subplot(312);
pcolor(RBR2.RBRgridtot.time,RBR2.RBRgridtot.z,RBR2.RBRgridtot.T);
shading flat 
axis ij
ylabel('depth')
xlabel('time')
datetick

ax(3)=subplot(313);
pcolor(RBRgrid.time,RBRgrid.z,RBRgrid.T);
shading flat 
axis ij
ylabel('depth')
xlabel('time')
datetick

linkaxes(ax)

figure
pcolor(RBRgrid.time,RBRgrid.z,RBRgrid.T);
shading flat 
axis ij
ylabel('depth','fontsize',15)
xlabel('time','fontsize',15)
set(gca,'fontsize',15)
datetick
fig=gcf;
fig.PaperPosition = [0 0 20 10];
print('-dpi','JohnWesleyPowell_RBR_T.png','-dpng2')

figure
pcolor(RBRgrid.time,RBRgrid.z,log10(RBRgrid.BScat));
shading flat 
axis ij
ylabel('depth','fontsize',15)
xlabel('time','fontsize',15)
set(gca,'fontsize',15)
datetick
title('log10(BScat)','fontsize',15)
print('JohnWesleyPowell_RBR_BScat.png','-dpng2')


