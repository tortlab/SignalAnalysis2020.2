% List of exercises 04

clear, clc, clf

load('LFP_HG_HFO.mat')

srate = 1000;
dt = 1/srate;
t = dt*(1:length(lfpHG));

%% Visualizing the raw signal

plot(t,lfpHFO)
hold on
plot(t,lfpHG)
hold off

%% Question 1

thetaA = eegfilt(lfpHG,srate,5,10);
thetaB = eegfilt(lfpHFO,srate,5,10);

gammaA = eegfilt(lfpHG,srate,60,100);
gammaB = eegfilt(lfpHFO,srate,60,100);

%% Question 2

AmpThetaA = abs(hilbert(thetaA));
AmpThetaB = abs(hilbert(thetaB));
AmpGammaA = abs(hilbert(gammaA));
AmpGammaB = abs(hilbert(gammaB));

%% Question 3

subplot(211)
plot(t,thetaA,'b-')
hold on
plot(t,AmpThetaA,'b-','linew',2)
plot(t,thetaB,'r-')
plot(t,AmpThetaB,'r-','linew',2)
hold off
xlim([1 2])

subplot(212)
plot(t,gammaA,'b-')
hold on
plot(t,AmpGammaA,'b-','linew',2)
plot(t,gammaB,'r-')
plot(t,AmpGammaB,'r-','linew',2)
hold off
xlim([1 2])

%% Question 4

subplot(221)
plot(AmpThetaA,AmpThetaB,'ko')
[r p] = corr(AmpThetaA',AmpThetaB')
title(['r = ' num2str(r)])

subplot(222)
plot(AmpThetaA,AmpGammaA,'ko')
[r p] = corr(AmpThetaA',AmpGammaA')
title(['r = ' num2str(r)])

subplot(223)
plot(AmpThetaB,AmpGammaB,'ko')
[r p] = corr(AmpGammaB',AmpThetaB')
title(['r = ' num2str(r)])

subplot(224)
plot(AmpGammaA,AmpGammaB,'ko')
[r p] = corr(AmpGammaA',AmpGammaB')
title(['r = ' num2str(r)])

%% Question 5

plot(t,lfpHG,'k-')
hold on
plot(t,thetaA-1,'b-')
plot(t,AmpThetaA-1,'b-','linew',2)
plot(t,gammaA-1.5,'r-')
plot(t,AmpGammaA-1.5,'r-','linew',2)

plot([40.5 40.7],[-1.8 -1.8],'k-','linew',2)
plot([40.5 40.5],[-1.8 -1.6],'k-','linew',2)

hold off
box off

xlim([40 44])
axis off

%% Question 6

freqvector = 1:1:100
bandwidth = 4

clear AmpSpectrum*
for flow = freqvector

fhigh = flow+bandwidth;
filtered = eegfilt(lfpHG,srate,flow,fhigh,0,500);
% filtered = eegfilt(LFP,srate,flow,fhigh);
AmpSpectrumHG(flow) = mean(abs(hilbert(filtered)));

filtered = eegfilt(lfpHFO,srate,flow,fhigh,0,500);
% filtered = eegfilt(LFP,srate,flow,fhigh);
AmpSpectrumHFO(flow) = mean(abs(hilbert(filtered)));
end

% %
clf
plot(freqvector+bandwidth/2,AmpSpectrumHG,'r-')
hold on
plot(freqvector+bandwidth/2,AmpSpectrumHFO,'b-')
hold off
xlabel('Frequency (Hz)')
ylabel('Amplitude (mV)')

%% Question 7

clear TFD*

freqvector = 0:1:100;

tic
count = 0;
for f = freqvector
    count = count + 1;
    f
    
   LFPfiltered = eegfilt(lfpHG,srate,f,f+4,0,300);
   AmpEnv = abs(hilbert(LFPfiltered));
   TFDHG(count,:) = AmpEnv;

      LFPfiltered = eegfilt(lfpHFO,srate,f,f+4,0,300);
   AmpEnv = abs(hilbert(LFPfiltered));
   TFDHFO(count,:) = AmpEnv;


end
toc

figure(2)

subplot(211)
imagesc(t,freqvector+2,TFDHG)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar
caxis([0 0.5])

subplot(212)
imagesc(t,freqvector+2,TFDHFO)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar

%% Question 8

Fc = 1;
NCycles = 7;
sd = NCycles/(2*pi*Fc);
Fb = 2*sd^2;
wname = ['cmor' num2str(Fb) '-' num2str(Fc)];
Freq = 1:0.5:100;
% scale = centfreq(wname)./(Freq*dt);
scale = Fc./(Freq*dt);

tic
WTHG = cwt(lfpHG,scale,wname);
WTHFO = cwt(lfpHFO,scale,wname);

toc

% %
TFDHG = abs(WTHG);
TFDHFO = abs(WTHFO);

subplot(211)
imagesc(t,Freq,TFDHG)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar

subplot(212)
imagesc(t,Freq,TFDHFO)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar




