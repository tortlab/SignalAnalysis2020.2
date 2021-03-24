% Lesson 06 - Aliasing and Nyquist frequency

clc, clear, clf

srate = 10000;
dt = 1/srate;
t = 0:dt:3;

% freq = 550;
% signal = sin(2*pi*freq*t+pi/2);

signal = sin(2*pi*120*t+pi/2)+sin(2*pi*720*t+pi/2);


subplot(211)
plot(t,signal)
xlim([0 0.05])

win = 1*srate;
nfft = 4*srate;
[Pxx F] = pwelch(signal,win,[],nfft,srate);

subplot(212)
plot(F,Pxx,'linew',3)
xlim([0 800])
xlabel('Freq (Hz)'),ylabel('Power')

% % sub-sampling

samplingfactor = 10;
signal2 = signal(1:samplingfactor:end);
srate2 = srate/samplingfactor;
t2 = t(1:samplingfactor:end);

subplot(211)
hold on
plot(t2,signal2,'ko-')
hold off

% %

win2 = 1*srate2;
nfft2 = 4*srate2;
[Pxx2 F] = pwelch(signal2,win2,[],nfft2,srate2);

subplot(212)
hold on
plot(F,Pxx2,'k-','linew',3)
hold off
xlim([0 1500])
xlabel('Freq (Hz)'),ylabel('Power')

%%

clc, clear, clf

srate = 10000;
dt = 1/srate;
t = 0:dt:3;

freq = 700;
signal = sin(2*pi*freq*t+pi/2);

% signal = sin(2*pi*120*t+pi/2)+sin(2*pi*720*t+pi/2);

subplot(211)
plot(t,signal)
xlim([0 0.05])

win = 1*srate;
nfft = 4*srate;
[Pxx F] = pwelch(signal,win,[],nfft,srate);

subplot(212)
plot(F,Pxx,'linew',3)
xlim([0 800])
xlabel('Freq (Hz)'),ylabel('Power')

% % sub-sampling

samplingfactor = 10;
signal2 = signal(1:samplingfactor:end);
srate2 = srate/samplingfactor;
t2 = t(1:samplingfactor:end);

subplot(211)
hold on
plot(t2,signal2,'ko-')
hold off
title(['Sampling rate = ' num2str(srate2) ' Hz'])

% %

win2 = 1*srate2;
nfft2 = 4*srate2;
[Pxx2 F] = pwelch(signal2,win2,[],nfft2,srate2);

subplot(212)
hold on
plot(F,Pxx2,'k-','linew',3)
hold off
xlim([0 1500])
xlabel('Freq (Hz)'),ylabel('Power')

%% Fixing aliasing effects

% pro tip: before subsampling, filter your signal
% below Nyquist frequency

FNyquist = srate2/2;

filtered = eegfilt(signal,srate,0,FNyquist);
subplot(211) 
hold on
plot(t,filtered,'r-')
hold off
xlim([0.01 0.06])


%%
[Pxx F] = pwelch(filtered,win,[],nfft,srate);

subplot(212)
hold on
plot(F,Pxx,'r-','linew',3)
xlim([0 1500])
xlabel('Freq (Hz)'),ylabel('Power')

% subsampling the filtered signal
samplingfactor = 10;
signal2 = filtered(1:samplingfactor:end);
srate2 = srate/samplingfactor;
t2 = t(1:samplingfactor:end);

subplot(211)
hold on
plot(t2,signal2,'ko--')
hold off

%%
[Pxx2 F] = pwelch(signal2,win2,[],nfft2,srate2);

subplot(212)
hold on
plot(F,Pxx2,'c-','linew',3)
hold off
xlim([0 1500])
xlabel('Freq (Hz)'),ylabel('Power')
