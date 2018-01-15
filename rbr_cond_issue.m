%load RBRprofiles
k=1;
clear T; clear C
for i=107:117
    if length(RBRprofiles{i}.T)<659
        keyboard
    end
    i
    T(:,k)=RBRprofiles{i}.T(10:489);
    C(:,k)=RBRprofiles{i}.C(10:489);
    k=k+1;
end

%%
dt=1/6; %cycles per second
clear specsT specsC
clear sd
clear fd
k=1;

for prf=1:11
    data=T(:,prf);
    seg_length=6*10; %20 sec worth of data
    N=length(data);
    M=seg_length/2;
    
    for n=1:floor(N/M-1)
        d=data((n-1)*M+1:(n+1)*M); %select data for the nth segment
        fd(:,n)=fft(d.*hanning(2*M)');
    end
    sd=sum(abs(fd(1:M+1,:)).^2,2)/N; %sum over all spectra (2nd index)
    sd(2:end)=sd(2:end)*2;
    
    df=1/(seg_length*dt);	%	cycles per dt
    freq=[0:M]*df;
    specsT(k,:)=sd';
    k=k+1;
end

for prf=1:11
    data=C(:,prf);
    seg_length=6*10; %20 sec worth of data
    N=length(data);
    M=seg_length/2;
    
    for n=1:floor(N/M-1)
        d=data((n-1)*M+1:(n+1)*M); %select data for the nth segment
        fd(:,n)=fft(d.*hanning(2*M)');
    end
    sd=sum(abs(fd(1:M+1,:)).^2,2)/N; %sum over all spectra (2nd index)
    sd(2:end)=sd(2:end)*2;
    
    df=1/(seg_length*dt);	%	cycles per dt
    freq=[0:M]*df;
    specsC(k,:)=sd';
    k=k+1;
end

%%
SN_figure(747474); clf;
h1(1)=loglog(freq(2:end),nanmean(specsC(:,2:end))./nanmean(specsC(:,3)),'.-')
hold on
h1(2)=loglog(freq(2:end),nanmean(specsT(:,2:end))./nanmean(specsT(:,3)),'.-')
set(h1,'linewidth',2)
xlabel('Freq (Hz)')
grid on
ylabel('Power Spectral Density ((mS/cm)^2 / °C^2 Hz^{-1})')
legend('Conductivity','Temperature')
title('Comparison of Conductivity and temperature spectra from Concerto 066046')

%%

Tprofiles{1}=RBRprofiles{104}.T
Tprofiles{2}=RBRprofiles{105}.T
Tprofiles{3}=RBRprofiles{107}.T
Tprofiles{4}=RBRprofiles{108}.T
Tprofiles{5}=RBRprofiles{109}.T
Tprofiles{6}=RBRprofiles{110}.T
Tprofiles{7}=RBRprofiles{111}.T
Tprofiles{8}=RBRprofiles{112}.T
Tprofiles{9}=RBRprofiles{113}.T
Tprofiles{10}=RBRprofiles{114}.T
Tprofiles{11}=RBRprofiles{115}.T
Tprofiles{12}=RBRprofiles{116}.T
Tprofiles{13}=RBRprofiles{117}.T
Tprofiles{14}=RBRprofiles{118}.T
Tprofiles{15}=RBRprofiles{119}.T

Cprofiles{1}=RBRprofiles{104}.C
Cprofiles{2}=RBRprofiles{105}.C
Cprofiles{3}=RBRprofiles{107}.C
Cprofiles{4}=RBRprofiles{108}.C
Cprofiles{5}=RBRprofiles{109}.C
Cprofiles{6}=RBRprofiles{110}.C
Cprofiles{7}=RBRprofiles{111}.C
Cprofiles{8}=RBRprofiles{112}.C
Cprofiles{9}=RBRprofiles{113}.C
Cprofiles{10}=RBRprofiles{114}.C
Cprofiles{11}=RBRprofiles{115}.C
Cprofiles{12}=RBRprofiles{116}.C
Cprofiles{13}=RBRprofiles{117}.C
Cprofiles{14}=RBRprofiles{118}.C
Cprofiles{15}=RBRprofiles{119}.C

Pprofiles{1}=RBRprofiles{104}.P
Pprofiles{2}=RBRprofiles{105}.P
Pprofiles{3}=RBRprofiles{107}.P
Pprofiles{4}=RBRprofiles{108}.P
Pprofiles{5}=RBRprofiles{109}.P
Pprofiles{6}=RBRprofiles{110}.P
Pprofiles{7}=RBRprofiles{111}.P
Pprofiles{8}=RBRprofiles{112}.P
Pprofiles{9}=RBRprofiles{113}.P
Pprofiles{10}=RBRprofiles{114}.P
Pprofiles{11}=RBRprofiles{115}.P
Pprofiles{12}=RBRprofiles{116}.P
Pprofiles{13}=RBRprofiles{117}.P
Pprofiles{14}=RBRprofiles{118}.P
Pprofiles{15}=RBRprofiles{119}.P

Sprofiles{1}=RBRprofiles{104}.S
Sprofiles{2}=RBRprofiles{105}.S
Sprofiles{3}=RBRprofiles{107}.S
Sprofiles{4}=RBRprofiles{108}.S
Sprofiles{5}=RBRprofiles{109}.S
Sprofiles{6}=RBRprofiles{110}.S
Sprofiles{7}=RBRprofiles{111}.S
Sprofiles{8}=RBRprofiles{112}.S
Sprofiles{9}=RBRprofiles{113}.S
Sprofiles{10}=RBRprofiles{114}.S
Sprofiles{11}=RBRprofiles{115}.S
Sprofiles{12}=RBRprofiles{116}.S
Sprofiles{13}=RBRprofiles{117}.S
Sprofiles{14}=RBRprofiles{118}.S
Sprofiles{15}=RBRprofiles{119}.S

%%

% filter T with a 1 Hz low pass


[b,a]=butter(4,0.5/3);

for i=1:15
    Tfilt{i}=filtfilt(b,a,Tprofiles{i})
    Cfilt{i}=filtfilt(b,a,Cprofiles{i})
end

for i=1:15
    Sclean{i}=gsw_SP_from_C(Cfilt{i},Tfilt{i},Pprofiles{i})
end

%%
for n=6;
n
SN_figure(2881); clf
plot(Sprofiles{n},Pprofiles{n})
hold on
plot(Sclean{n},Pprofiles{n})
axis ij
ylabel('Depth')
xlabel('Salinity')

SN_figure(2882); clf
plot(Tprofiles{n},Pprofiles{n})
hold on
plot(Tfilt{n},Pprofiles{n})
axis ij
ylabel('Depth')
xlabel('Temperature')


end












