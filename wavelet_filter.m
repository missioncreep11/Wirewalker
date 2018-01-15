function wavelet_filter(WWmeta)

    load([WWmeta.aqdpath 'Profiles_' WWmeta.name_aqd],'AQDprofiles')
    N=cellfun(@(x) length(x.Burst_MatlabTimeStamp),AQDprofiles,'un',0); % get length of all upcast
    %indOK=find([N{:}]>max([N{:}])/2);
    data=AQDprofiles(find([N{:}]>max([N{:}])/2));
%%anonymous
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();

%% cells of upcasts 
    Velocell=cellfun(@(x) complex(medfilt1(x.Burst_VelEast,20),...
                                  medfilt1(x.Burst_VelNorth,20)),...
                     data,'un',0);

%time=arrayfun(@(x) tu(datastruct.idx(x,1):datastruct.idx(x,2)),1:length(datastruct.idx),'un',0);


%% some parameters 
    w0   = 6.; % for the wavelet
    dt=mean(diff(data{1}.Burst_MatlabTimeStamp))*86400;
    Fs=1/dt;

    N=cellfun(@length,Velocell,'un',0); % get length of all upcast
    mu=cellfun(@mean,Velocell,'un',0);  % get the mean of all upcast
    Velop=cellfun(@(x,y) x-y,Velocell,mu,'un',0); % remove the mean of all upcast

    % compute wavelet decomposition for the real and imaginary part of the
    % velocities. I have to do it in 2 times (real and imaginary) since my
    % wavelet decomposition applied only on the positive part of the fourrier
    % transform
    [rcwt,~,~,~]=cellfun(@(x) morlet_wavelet(real(x),Fs,w0,[]),Velop,'un',0);
    [icwt,scales,conversion,coi]=cellfun(@(x) morlet_wavelet(imag(x),Fs,w0,[]),Velop,'un',0);

    
    % reconstruct the signal. I tuned a little bit with my .62 coef to make it
    % close to the original signal there is most likely a proper way to do it
    % check D. rudnick wavelet routines
    Cs=0.7784;r0=0.7511;%from theory
    ds=cellfun(@(x) nanmean(diff(x)),scales,'un',0);
    recons_coef=cellfun(@(x,y) (x*sqrt(dt)/Cs/r0./sqrt(y)),ds,scales,'un',0);
    r_recons_signal=cellfun(@(x,y,z) (.62*repmat(x,[1,z]).*abs(y).*cos(angle(y))),recons_coef,rcwt,N,'un',0);
    i_recons_signal=cellfun(@(x,y,z) (.62*repmat(x,[1,z]).*abs(y).*cos(angle(y))),recons_coef,icwt,N,'un',0);


    maxrecons=cellfun(@(x,y) max(mean(x,1)+real(y)),r_recons_signal,mu);
    maxsignal=cellfun(@(x,y) max(real(x(1:y))),Velocell,N);
    adjustcoef=maxsignal./maxrecons;

    if 1==0
        i=10
        figure(1)
        subplot(212)
        imagesc(datenum2yday(data{i}.Burst_MatlabTimeStamp),...
            scales{i},abs(rcwt{i}));
        caxis([0,1])
        colorbar;
        title('wavelet signal')
        ylabel('scales (TF*period)')
        xlabel('year day')
        
        
        subplot(211)
        hold on
        plot(datenum2yday(data{i}.Burst_MatlabTimeStamp),real(Velocell{i}(1:N{i})),'b')
        plot(datenum2yday(data{i}.Burst_MatlabTimeStamp),adjustcoef(i)*mean(r_recons_signal{i},1)+real(mu{i}),'r')
        title('wavelet  reconstruction signal ')
        legend('Raw data','total wavelet reconstruction',...
            'location','best')
        xlabel('time axis')
        ylabel('m/s')
        hold off
        
        figure(2)
        cmap=colormap(jet(200));
        ax(1)=subplot('Position',[0.1 0.5 0.5 .4]);
        %imagesc(abs(cwt));caxis([0,1]);colorbar
        imagesc(datenum2yday(data{i}.Burst_MatlabTimeStamp),scales{i},r_recons_signal{i});caxis([-1,1]);
        title('wavelet signal')
        ylabel('scales (i.e: Freq)')
        xlabel('time')
        ax(2)=subplot('Position',[0.1 0.1 0.5 .2]);
        plot(datenum2yday(data{i}.Burst_MatlabTimeStamp),std(r_recons_signal{i},[],1))
        xlabel('time')
        title('Scale std')
        xlim(datenum2yday(data{i}.Burst_MatlabTimeStamp([1,N{i}])))
        %plot(std(abs(cwt),[],1))
        ax(3)=subplot('Position',[0.7 0.5 0.2 .4]);
        plot(std(r_recons_signal{i},[],2),scales{i});axis ij
        %plot(std(abs(cwt),[],2),1:length(scales));axis ij
        ylim(scales{i}([1,length(scales{i})]))
        ylabel('Scales (TF*Period)')
        title('time std')
    end
    
    fc=1/30; %3 cpd ~ 3h
    fnb=1/(2*dt)  ;        % nyquist freq
    [b, a]   = butter(3,fc/fnb,'low');
    Fstdsig=cellfun(@(x) filtfilt(b,a,std(x,[],1)),r_recons_signal,'un',0);
    ddFstdsig=cellfun(@(x) diff(x,1,2),Fstdsig,'un',0);


    hold(ax(2),'on')
    plot(ax(2),datenum2yday(data{i}.Burst_MatlabTimeStamp),Fstdsig{i},'r','linewidth',2)
    hold(ax(2),'off')

    %% find the "depth" of the maximun oscilation (via the peak of the second
    %  derivative of the time standard deviation of the reconstructed signal)
    [M,I]=cellfun(@max, ddFstdsig,'un',0);
    % if no peak depth is nan
    D_surf_noise=cellfun(@(x,y,z) iif(isempty(find(y(1:z-1)<0,1,'last')),NaN,...
        ~isempty(find(y(1:z-1)<0,1,'last')),...
        x-find(y(1:z-1)<0,1,'last')),...
        N,ddFstdsig,I,'un',0 );

    % transform Depth into a vector
    D_surf_noise=cell2mat(D_surf_noise);
    % interp to remove nan
    ind_upcast=1:length(D_surf_noise);
    D_surf_noise=interp1(ind_upcast(~isnan(D_surf_noise)),D_surf_noise(~isnan(D_surf_noise)),...
        ind_upcast);
    if(sum(isnan(D_surf_noise))>0)
        D_surf_noise(1)=nanmean(D_surf_noise);
        D_surf_noise(end)=nanmean(D_surf_noise);
    end
    if 1==0
        hold(ax(2),'on')
        plot(ax(2),datenum2yday(data{i}.Burst_MatlabTimeStamp(end-D_surf_noise([i i]))),[min(Fstdsig{i}) max(Fstdsig{i})],'k','linewidth',2)
        hold(ax(2),'off')
        figure(3)
        plot(1:length(D_surf_noise),D_surf_noise)
        xlabel('nb upcast')
        ylabel('Depth of signal wave (in index)')
    end

    % design a filter which gonna smooth the Depth vector every 10 upcast
    fc=1/10; %cut frequency every 10 upcast
    fnb=1/(2*dt)  ;        % nyquist freq
    [b, a]   = butter(3,fc/fnb,'low');
    D_surf_noise=filtfilt(b,a,D_surf_noise);
    
    % arbitrary increase the depth to be sure that we remove the suface wave
    % signal deeply enough
    D_surf_noise=ceil(2*D_surf_noise); % watch possible cases where 2*D_surf_noise longer than N (not good); to fix
    N=cell2mat(N);
    D_surf_noise(D_surf_noise>N)=N(D_surf_noise>N)-2;
    N=num2cell(N);
    D_surf_noise=num2cell(D_surf_noise);
    if 1==0
        hold on
        plot(1:length(D_surf_noise),[D_surf_noise{:}],'r','linewidth',2)
        legend('real surface signal depth','2*smooth depth signal')
        hold off
    end
        
        
    
    
    %% find the width of the frequency band of the surface wave signal using
    %  (via the peak of the second derivative of the scale standard deviation
    %  of the reconstructed signal)
    ddstd2sig=cellfun(@(x,y) diff(std(x(:,end-ceil(y):end),[],2)',1,2),...
        r_recons_signal,D_surf_noise,'un',0);
    
    % scale 1:N/2 pretty safe place to look at surface wave signal
    % would make more sens if I transform scale to frequency. scale 27 ~ freq >
    % few minutes
    [M,I]=cellfun(@(x) max(std(x(1:floor(size(x,1)/2),:),1,2)),r_recons_signal,'un',0);
    
    
    S_surf_noise=cellfun(@(x,y) iif(isempty(find(x(1:y-1)<0,1,'last')),NaN,...
        ~isempty(find(x(1:y-1)<0,1,'last')),...
        find(x(1:y-1)<0,1,'last')),...
        ddstd2sig,I,'un',0);
    
    S_surf_noise=cell2mat(S_surf_noise);
    ind_upcast=1:length(S_surf_noise);
    S_surf_noise=interp1(ind_upcast(~isnan(S_surf_noise)),S_surf_noise(~isnan(S_surf_noise)),...
        ind_upcast);
    I=cell2mat(I);
    
    if ~isempty(isnan(S_surf_noise))
        S_surf_noise(isnan(S_surf_noise))=...
            I(isnan(S_surf_noise))-floor(nanmean(I-S_surf_noise));
    end
    if any(I==1)
        S_surf_noise(I==1)=floor(mean(S_surf_noise));
        I(I==1)=ceil(mean(I));
    end


    
    % design a filter which gonna smooth the Depth vector every 10 upcast
    fc=1/30; %cut frequency every 10 upcast
    dtup=7*60; % one upcast every 7 min (approximation)
    fnb=1/(2*dt)  ;        % nyquist freq
    [b, a]   = butter(3,fc/fnb,'low');
    S_surf_noise=ceil(filtfilt(b,a,S_surf_noise));
    I_surf_noise=ceil(filtfilt(b,a,I));
    I_surf_noise=num2cell(I_surf_noise);
    S_surf_noise=num2cell(S_surf_noise);
    
    W_surf_noise=cellfun(@(x,y) (x-y),I_surf_noise,S_surf_noise,'un',0);
    gaussfilt=cellfun(@(x,y) (x*gausswin(ceil(2*y))),D_surf_noise,W_surf_noise,'un',0);
    scalefilt=cellfun(@(x) zeros(length(x),1),scales,'un',0);
    ffilt1=struct([]);ffilt2=struct([]);ffilt3=struct([]);
    r_ffilt1cwt=struct([]);r_ffilt2cwt=struct([]);r_ffilt3cwt=struct([]);
    i_ffilt1cwt=struct([]);i_ffilt2cwt=struct([]);i_ffilt3cwt=struct([]);
    r_ffiltcwt=struct([]);
    i_ffiltcwt=struct([]);
    r_recons_filtsignal=struct([]);
    i_recons_filtsignal=struct([]);
for i=1:length(S_surf_noise)
    %disp(i)
    scalefilt{i}(S_surf_noise{i}:S_surf_noise{i}+ceil(2*W_surf_noise{i}-1))=gaussfilt{i};
    time_ind=1:N{i};
    [ttind,gauss2D]=meshgrid(time_ind,scalefilt{i});
    filterspace=min(length(scales{i}),1.3*W_surf_noise{i})/N{i}*time_ind+S_surf_noise{i}+min(length(scales{i}),1.3*W_surf_noise{i});
    [filterspace2D,sscale]=meshgrid(filterspace,1:length(scales{i}));
    
    ffilt1{i}=ttind>=gauss2D;
    ffilt1{i}=flipud(ffilt1{i}')';
    ffilt2{i}=sscale<=filterspace2D;
    ffilt3{i}=ffilt1{i};
    ffilt3{i}(ffilt2{i}==0)=0;

    r_ffilt1cwt{i}=ffilt1{i}.*rcwt{i};
    r_ffilt2cwt{i}=~ffilt2{i}.*rcwt{i};
    r_ffilt3cwt{i}=ffilt3{i}.*rcwt{i};
    i_ffilt1cwt{i}=ffilt1{i}.*icwt{i};
    i_ffilt2cwt{i}=~ffilt2{i}.*icwt{i};
    i_ffilt3cwt{i}=ffilt3{i}.*icwt{i};

    
    r_background_max=max(abs(r_ffilt3cwt{i}(:)));
    i_background_max=max(abs(i_ffilt3cwt{i}(:)));
    r_background_level=.33*r_background_max;
    i_background_level=.33*i_background_max;
    
    r_ffilt3cwt{i}=r_background_level*r_ffilt3cwt{i}./r_background_max;
    i_ffilt3cwt{i}=i_background_level*i_ffilt3cwt{i}./i_background_max;
    
    r_ffiltcwt{i}=r_ffilt2cwt{i};
    r_ffiltcwt{i}(abs(r_ffilt3cwt{i})~=0)=r_ffilt3cwt{i}(abs(r_ffilt3cwt{i})~=0);
    i_ffiltcwt{i}=i_ffilt2cwt{i};
    i_ffiltcwt{i}(abs(i_ffilt3cwt{i})~=0)=i_ffilt3cwt{i}(abs(i_ffilt3cwt{i})~=0);

    r_recons_filtsignal{i}=.62*repmat(recons_coef{i},[1,N{i}]).*abs(r_ffiltcwt{i}).*cos(angle(r_ffiltcwt{i}));
    i_recons_filtsignal{i}=.62*repmat(recons_coef{i},[1,N{i}]).*abs(i_ffiltcwt{i}).*cos(angle(i_ffiltcwt{i}));
end

for i=11:20
    figure(i)
    data1=abs(r_ffiltcwt{i});
    data1(data1==0)=nan;
    subplot('Position',[0.1 .55 .7 .4])
    pcolor(datenum2yday(data{i}.Burst_MatlabTimeStamp),1:length(scales{i}),data1);caxis([0,1]);
    shading interp
    axis ij 
    title('wavelet decomposition and automated nan filter')
    xlabel('year day')
    ylabel('Wavelet Scale index')
    cbar=colorbar;
    cbar.Position=[0.85 .55 .02 .4];
    
    subplot('Position',[0.1 .1 .7 .35])
    hold on
    plot(datenum2yday(data{i}.Burst_MatlabTimeStamp),real(Velocell{i}(1:N{i})),'b')
    plot(datenum2yday(data{i}.Burst_MatlabTimeStamp),adjustcoef(i).*...
                                (mean(r_recons_signal{i},1))+real(mu{i}),'r')
    plot(datenum2yday(data{i}.Burst_MatlabTimeStamp),adjustcoef(i).*...
                               (nanmean(r_recons_filtsignal{i},1))+real(mu{i}),'g')
    hold off
    title('wavelet  reconstruction signal ')
    legend('Raw data','total wavelet reconstruction','wavelet filter reconstruction',...
            'location','best')
    xlabel('time axis')
    xlim(datenum2yday(data{i}.Burst_MatlabTimeStamp([1,N{i}])))
    ylabel('m/s')
    %saveas(gcf,sprintf('%sWAVELET_FILTER/wavelet_filter_upcast%i_v2.png',figure_path,i))
    pause
end


adjustcoef=num2cell(adjustcoef);

velfilt=cellfun(@(x,y,z,xx) xx.*complex(nanmean(x,1)+real(z),nanmean(y,1)+imag(z)),...
                        r_recons_filtsignal,i_recons_filtsignal,mu,adjustcoef,'un',0);

subplot(211)
plot(eu,'k')
hold on
plot(real([velfilt{:}]),'r')
subplot(212)
plot(nu,'k')
hold on
plot(imag([velfilt{:}]),'r')
                    
% grid the filtered velocities
[m,n]=size(datastruct.std_profiles.T);
zgrid=0:.25:200;
datastruct.e_aqd_filt=real([velfilt{:}]);
datastruct.n_aqd_filt=imag([velfilt{:}]);

datastruct.std_profiles_filt=datastruct.std_profiles;
datastruct.std_profiles_filt.aqd=make_standard_profiles_AQDII_filt_surf(datastruct,datastruct.idx,zgrid); %110 max pressure d1

datastruct.std_profiles_filt.e_abs=datastruct.std_profiles_filt.aqd.e+repmat(datastruct.std_profiles.e_drift,m,1);
datastruct.std_profiles_filt.n_abs=datastruct.std_profiles_filt.aqd.n+repmat(datastruct.std_profiles.n_drift,m,1);
datastruct.std_profiles_filt.u_abs=datastruct.std_profiles_filt.aqd.u;

aquarius_d4=datastruct;
save([ROOT 'aquarius/WW/d4/rbr/aquarius_d4.mat'],'aquarius_d4')

