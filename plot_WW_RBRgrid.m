function [hf,ha] = plot_WW_RBRgrid(gridin,begtime,endtime,vars,hf);
% Function to make nice plots of gridded WW data.
% possible values for vars:
% From RBR: 'T','S','C','BScat','F_chla','F_CDOM','DO'
% From AQD: 'u','v','w','amp'
% From post-processing: 'n2','rho','L_ot','eps_ot','Lmin'
%
% in order to access all these variables, gridin should be EITHER:
% 1) WWgrid generated from get_WW_data, OR
% 2) WWmeta with the fields generated in process_PROJECTNAME.m
% 
% begtime and endtime are the range of days to plot (matlab datenum format)
%
% Created by M. Hamann 9/11/17

if nargin<4
    vars = {'T','S','BScat','F_chla','DO'};
end
if nargin<3 || isempty(endtime)
    endtime = gridin.time(end);
end
if nargin<2 || isempty(begtime)
    begtime = gridin.time(1);
end
if ~isfield(gridin,'time')
    gridin = get_WW_data(gridin,begtime,endtime);
end

id = find(gridin.time>=begtime & gridin.time<=endtime);

if size(gridin.z,2)~=1
    gridin.z = gridin.z';
end

isu = strcmp(vars,'u');
isu2 = find(isu == 1);
isv = strcmp(vars,'v');
isw = strcmp(vars,'w');
if ~isempty(isu2);
    isu3 = find(~isnan(gridin.u(:,id)));
    if isempty(isu3)
        vars2 = vars(isu==0 & isv==0 & isw==0);
        vars = vars2;
    end
end
n = length(vars);

if nargin<5
    hf = figure(153); clf; tallfigure;
else
    figure(hf);
end
ha = MySubplot(0.05,0.15,0,0.05,0.1,0.02,1,n);

for ii = 1:n
    axes(ha(ii));
    if strcmp(vars{ii},'eps_ot') || strcmp(vars{ii},'n2')
        pcolor(gridin.time(id),gridin.z,real(log10(gridin.(vars{ii})(:,id))));
    else
        pcolor(gridin.time(id),gridin.z,gridin.(vars{ii})(:,id));
    end
    shading flat; axis ij
    
    c = cbstay;
    
    if strcmp(vars{ii},'T'); caxis([8 20]); ylabel(c,'Temp, {}^\circ{C}'); colormap(gca,'jet')
    elseif strcmp(vars{ii},'S'); caxis([33 34]); ylabel(c,'Sal, psu'); colormap(gca,redblue)
    elseif strcmp(vars{ii},'rho'); caxis([1024 1026.5]); ylabel(c,'\rho, kg m^{-3}'); 
    elseif strcmp(vars{ii},'BScat'); caxis([100 500]); ylabel(c,'Backscatter'); colormap(gca,hot)
    elseif strcmp(vars{ii},'F_chla'); caxis([50 300]); ylabel(c,'Chl, volts');colormap(gca,copper)
    elseif strcmp(vars{ii},'DO'); caxis([40 100]); ylabel(c,'DO, %'); colormap(gca,'jet')
    elseif strcmp(vars{ii},'u'); caxis([-0.25 0.25]); ylabel(c,'u, m/s'); colormap(gca,redblue)
    elseif strcmp(vars{ii},'v'); caxis([-0.25 0.25]); ylabel(c,'v, m/s'); colormap(gca,redblue)
    elseif strcmp(vars{ii},'w'); caxis([-0.05 0.05]); ylabel(c,'w, m/s'); colormap(gca,redblue)
    elseif strcmp(vars{ii},'eps_ot'); caxis([-9 -5]); ylabel(c,'\epsilon_{OT}, W/kg'); colormap(gca,parula)
    elseif strcmp(vars{ii},'n2'); caxis([-5 -3]); ylabel(c,'N^2, s^{-2}'); colormap(gca,parula)
    else ylabel(c,vars{ii})
    end
    
    if ii == n
%         datetick
    else
        set(gca,'xticklabel',[])
    end
    if ii == 1
        title([datestr(begtime),' to ',datestr(endtime)])
    end
end


linkaxes(ha,'x')
axes(ha(n)); datetick