% aula 10

% Hilbert transform

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:3;

LFP = sin(2*pi*8*t);

whitenoise = randn(size(t));
LFP = eegfilt(whitenoise,srate,8,10);

% Matlab's built-in 'hilbert' function outputs
% the analytical representation of a signal
% defined by: AR = LFP + iH(LFP), where
% H(LFP) is the Hilbert transform of the signal

AR = hilbert(LFP);

original_signal = real(AR);
HT = imag(AR);


%%

subplot(1,1,1)

plot(t,LFP)
hold on
% plot(t,original_signal)
plot(t,HT,'g-')
hold off


%% Amplitude envelope

Amp = abs(hilbert(LFP));

for n = 1:3000 % samples
    subplot(311)
    plot(t,LFP)
    hold on
    plot(t,HT,'r-')
    plot(t,Amp,'k-','linewidth',2)
    
    plot(t(n),LFP(n),'bo')
    plot(t(n),HT(n),'ro')
    plot(t(n),Amp(n),'ko')
    
    hold off
    
    
    xlim([0 3])
    
    subplot(3,1,[2 3])
    
    plot([0 LFP(n)],[0 HT(n)],'k-')
    hold on
    plot(LFP(n),HT(n),'ko')
    
    title(['Instantaneous amplitude = ' ...
        num2str(Amp(n))])
    
    axis square
    
    xlim([-.2 .2])
    ylim(xlim())
    
    
    hold off
    pause(0.01)
end


%%

LFP = sin(2*pi*8*t).*sin(2*pi*0.5*t);

HT = imag(hilbert(LFP));

AmpEnv = abs(hilbert(LFP));

subplot(311)
plot(t,LFP)
hold on
plot(t,AmpEnv,'k-','linew',3)
hold off

subplot(3,1,[2 3])
% v = 1:10;
% plot3(v,v.^2,v) % plot3(x,y,z)

plot3(t,LFP,HT)
hold on
plot3(t,zeros(size(t)),zeros(size(t)),...
    'r-','linew',1)
hold off

xlabel('time (s)')
ylabel('Real')
zlabel('Imag')

%% obs: 

% the hilbert transform of a constant signal is zero

X = ones(1,1000);
HT = imag(hilbert(X));

AmplitudeEnvelope = abs(hilbert(X))



%% Amplitude Spectrum

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:10;

LFP = 2*sin(2*pi*20*t) + sin(2*pi*70*t);

order = 300 % Filter order

f = 10;

kernel = sin(2*pi*f*t(1:order));
kernel = kernel/sum(kernel.^2);

Filtered = conv(LFP,kernel,'same');
Filtered = conv(Filtered,kernel(end:-1:1),'same');

HT = imag(hilbert(Filtered));
Amp = abs(hilbert(Filtered));

% mean amplitude
MeanAmp = mean(Amp);

% Root Mean Square
RMS = sqrt(mean(Filtered.^2));

subplot(311)
plot(t,LFP)
hold on
plot(t,Filtered,'b-')
plot(t,HT,'r-')
plot(t,Amp,'k-','linew',2)
xlim([1 1.2])
hold off

title(['Freq = ' int2str(f) ...
    ' Hz;  Amp = ' num2str(MeanAmp) ...
    ';  RMS = ' num2str(RMS)])

%% Iterating over frequencies

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:10;

% LFP = 2*sin(2*pi*20*t) + 1*sin(2*pi*70*t);

LFP = sin(2*pi*20*t) + sin(2*pi*40*t)+ sin(2*pi*70*t);

order = 200 % Filter order
freqvector = 1:1:100;

for f = freqvector
    
    kernel = sin(2*pi*f*t(1:order));
    kernel = kernel/sum(kernel.^2);
    
    Filtered = conv(LFP,kernel,'same');
    Filtered = conv(Filtered,kernel(end:-1:1),'same');
    
    HT = imag(hilbert(Filtered));
    Amp = abs(hilbert(Filtered));
    
    % mean amplitude
    MeanAmp = mean(Amp);
    AmpSpectrum(f) = MeanAmp;
    
    % Root Mean Square
    RMS = sqrt(mean(Filtered.^2));
    RMSSpectrum(f) = RMS;
    
    
    subplot(311)
    plot(t,LFP)
    hold on
    plot(t,Filtered,'b-')
    plot(t,HT,'r-')
    plot(t,Amp,'k-','linew',2)
    xlim([1 1.2])
    hold off
    
    title(['Freq = ' int2str(f) ...
        ' Hz;  Amp = ' num2str(MeanAmp) ...
        ';  RMS = ' num2str(RMS)])
    
    pause(0.05)
    
end

% %

subplot(3,2,[3 5])
plot(freqvector,AmpSpectrum)
xlabel('Frequency (Hz)')
ylabel('Amplitude (mV)')

subplot(3,2,[4 6])
plot(freqvector,RMSSpectrum)
xlabel('Frequency (Hz)')
ylabel('RMS (mV)')


%% Amplitude spectrum via eegfilt

srate = 1000;
dt = 1/srate;
t = dt:dt:10;

LFP = sin(2*pi*20*t) + sin(2*pi*40*t)+ sin(2*pi*70*t);

freqvector = 1:1:100
bandwidth = 4
clear AmpSpectrum

for flow = freqvector
    
    fhigh = flow+bandwidth;
    filtered = eegfilt(LFP,srate,flow,fhigh);
    AmpSpectrum(flow) = mean(abs(hilbert(filtered)));
    
end

% %
% subplot(3,2,[4 6])
% subplot(111)
subplot(211)
plot(freqvector+bandwidth/2,AmpSpectrum)
xlabel('Frequency (Hz)')
ylabel('Amplitude (mV)')

subplot(212)
[Pxx F] = pwelch(LFP,2*srate,[],freqvector,srate)
plot(F,Pxx)
xlim([0 120])


%% TFD - time frequency decomposition

LFP = sin(2*pi*10*t) + sin(2*pi*40*t);
LFP(1:5000) = 0;

% LFPfiltered  = eegfilt(LFP,srate,6.3,15.87234);
% AmpEnv = abs(hilbert(LFPfiltered));

% plot(t,LFP)
% hold on
% plot(t,LFPfiltered,'r-')
% plot(t,AmpEnv,'k-','linew',2)
% hold off
% 
% xlim([4.5 5.5])

% %

clear TFD

freqvector = 0:1:100;

tic
count = 0;
for f = freqvector
    count = count + 1;
    
   LFPfiltered = eegfilt(LFP,srate,f,f+4,0,300);
   AmpEnv = abs(hilbert(LFPfiltered));
   TFD(count,:) = AmpEnv;
   
% %   plot(t,LFP)
% % hold on
% % plot(t,LFPfiltered,'r-')
% % plot(t,AmpEnv,'k-','linew',2)
% % hold off
% % 
% % title(['Frequency (Hz) = ' num2str(f+2) ' Hz'])
% % xlim([4.5 5.5])
% % 
% % pause(0.01)
   
end
toc


%% Plotting TFD

subplot(211)
imagesc(t,freqvector+2,TFD)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar


subplot(212)
tic
[S,F,T,TFD2]=spectrogram(LFP,2*srate,1.9*srate,freqvector+2,srate);
toc
imagesc(t,freqvector+2,TFD2)
axis xy
xlabel('Time (s)')
ylabel('Frequency (Hz)')
colorbar








