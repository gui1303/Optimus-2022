%% 1. Initialization 
clear all;close all;clc;

%% 2. Load Data
FileName                = 'IEA-15-240-RWT-UMaineSemiNoDamper_FB.outb';
FBnPD                   = ReadFASTbinaryIntoStruct(FileName);
FileName                = 'IEA-15-240-RWT-UMaineSemiWithDamper_FB.outb';
FBwPD                   = ReadFASTbinaryIntoStruct(FileName);

%% 3. Estimation of Spectrum from Data
%Platform pitch
signal               	= detrend(FBnPD.PtfmPitch,'constant');
nBlocks                 = 1;
nOverlap                = [];   % default: nDataPerBlock/2;
nFFT                    = [];   % default: 2^nextpow2(nDataPerBlock);
SamplingFrequency       = 1/diff(FBwPD.Time(1:2));
nDataPerBlock           = floor(size(signal,1)/nBlocks/2)*2; % should be even
vWindow                 = hamming(nDataPerBlock);
[S_FBnPD_PtfmPitch,f]  	= pwelch(signal,vWindow,nOverlap,nFFT,SamplingFrequency); 
signal         = detrend(FBwPD.PtfmPitch,'constant');
[S_FBwPD_PtfmPitch,f]  	= pwelch(signal,vWindow,nOverlap,nFFT,SamplingFrequency); 

% Rotor speed
signal                  = detrend(FBnPD.RotSpeed,'constant');
nBlocks                 = 1;
nOverlap                = [];   % default: nDataPerBlock/2;
nFFT                    = [];   % default: 2^nextpow2(nDataPerBlock);
SamplingFrequency       = 1/diff(FBwPD.Time(1:2));
nDataPerBlock           = floor(size(signal,1)/nBlocks/2)*2; % should be even
vWindow                 = hamming(nDataPerBlock);
[S_FBnPD_RotSpeed,f]  	= pwelch(signal,vWindow,nOverlap,nFFT,SamplingFrequency); 
signal                  = detrend(FBwPD.RotSpeed,'constant');
[S_FBwPD_RotSpeed,f]  	= pwelch(signal,vWindow,nOverlap,nFFT,SamplingFrequency); 

% Blade pitch
signal                  = detrend(FBnPD.BldPitch1,'constant');
nBlocks                 = 1;
nOverlap                = [];   % default: nDataPerBlock/2;
nFFT                    = [];   % default: 2^nextpow2(nDataPerBlock);
SamplingFrequency       = 1/diff(FBwPD.Time(1:2));
nDataPerBlock           = floor(size(signal,1)/nBlocks/2)*2; % should be even
vWindow                 = hamming(nDataPerBlock);
[S_FBnPD_BldPitch,f]  	= pwelch(signal,vWindow,nOverlap,nFFT,SamplingFrequency); 
signal                  = detrend(FBwPD.BldPitch1,'constant');
[S_FBwPD_BldPitch,f]  	= pwelch(signal,vWindow,nOverlap,nFFT,SamplingFrequency); 

%% 4. Plots
% frequency
figure('Name','Spectra')
hold all; grid on; box on
plot(f,S_FBnPD_PtfmPitch)
plot(f,S_FBwPD_PtfmPitch)
set(gca,'xScale','log')
set(gca,'yScale','log')
xlim([1e-3 1e0])
xlabel('frequency [Hz]')
ylabel('spectra [(deg)^2/Hz]')
legend('no Platform Damper','with Platform Damper')
title('Spectrum Platform Pitch Angle')

figure('Name','Spectra')
hold all; grid on; box on
plot(f,S_FBnPD_RotSpeed)
plot(f,S_FBwPD_RotSpeed)
set(gca,'xScale','log')
set(gca,'yScale','log')
xlim([1e-3 1e0])
xlabel('frequency [Hz]')
ylabel('spectra [(deg)^2/Hz]')
legend('no Platform Damper','with Platform Damper')
title('Spectrum Rotor Speed')

figure('Name','Spectra')
hold all; grid on; box on
plot(f,S_FBnPD_BldPitch)
plot(f,S_FBwPD_BldPitch)
set(gca,'xScale','log')
set(gca,'yScale','log')
xlim([1e-3 1e0])
xlabel('frequency [Hz]')
ylabel('spectra [(deg)^2/Hz]')
legend('no Platform Damper','with Platform Damper')
title('Blade Pitch Angle')
%%
% time
figure('Name','Platform Pitch Angle')
hold all; grid on; box on
plot(FBnPD.Time,FBnPD.PtfmPitch)
plot(FBwPD.Time,FBwPD.PtfmPitch)
xlabel('time [s]')
ylabel('PtfmPitch [deg]')


