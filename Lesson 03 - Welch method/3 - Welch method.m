%% Lesson 3 - Welch Method

% Reviewing lesson 2

srate = 1000; % in Hz
dt = 1/srate; % in s
f = 50; % in Hz
Tmax = 1;
t = dt:dt:Tmax;

X = sin(2*pi*f*t)+sin(2*pi*15*t);

subplot(3,1,1)
plot(t,X)

xlabel('Time (s)')
ylabel('mV')

subplot(3,1,[2 3])

ff = 50; % freq Fourier Kernell

K = exp(-1i*2*pi*ff*t);

plot(X.*K,'yo')

axis square

FX = mean(X.*K);

hold on
plot([0 real(FX)],[0 imag(FX)],'r-','linew',3)
hold off

xlim([-1 1])
ylim(xlim())

Power = FX*conj(FX);

title(['Power at ' num2str(ff) ' Hz = ' num2str(Power)])

%% Power Spectrum

srate = 100; % in Hz
dt = 1/srate; % in s
f = 15; % in Hz
Tmax = 1;
t = dt:dt:Tmax;

X = sin(2*pi*f*t)+0*sin(2*pi*15*t);

freq = 0:0.01:80;
clear PSD

count = 0;
for ff = freq % freq Fourier Kernell
count = count+1;
K = exp(-1i*2*pi*ff*t);
FX = sum(X.*K);
PSD(count) = FX*conj(FX)/Tmax;
end

figure(2)
subplot(1,1,1)
plot(freq,PSD,'k-')


hold on
I = find(freq==15.1);
plot(freq(I),PSD(I),'ro','markerf','m')
hold off

xlabel('Frequency (Hz)')
ylabel('Power')

xlim([10 20])

dF = freq(2)-freq(1)

%% Hamming window and Hann window

clf

srate = 1000; % in Hz
dt = 1/srate; % in s
f = 15; % in Hz
Tmax = 1;
t = dt:dt:Tmax;

X = sin(2*pi*f*t);

subplot(3,1,1)
plot(t,X)

subplot(3,1,2)
W1 = hann(length(X))';

% W1 = rectwin(length(X))';

plot(t,W1)

W2 = hamming(length(X))';
hold on
plot(t,W2)
hold off

subplot(3,1,3)
plot(t,X.*W1)
hold on
plot(t,X.*W2)
hold off

%% Hann and Hamming correction factors

freq = 0:0.01:80;
clear PSD PSDW1 PSDW2

count = 0;
for ff = freq % freq Fourier Kernell
count = count+1;
K = exp(-1i*2*pi*ff*t);
FX = sum(X.*K)*dt;
PSD(count) = FX*conj(FX)/Tmax;

FX = sum((X.*W1).*K)*dt*1.63; %length(W1)/sum(W1);
PSDW1(count) = FX*conj(FX)/Tmax;

FX = sum((X.*W2).*K)*dt*1.59; %length(W2)/sum(W2);
PSDW2(count) = FX*conj(FX)/Tmax;

end

figure(2)
subplot(1,1,1)
plot(freq,PSD,'k-')

hold on
plot(freq,PSDW1,'b-')
plot(freq,PSDW2,'r-')
hold off

hold on
I = find(freq==13.5);
plot(freq(I),PSD(I),'ro','markerf','m')
hold off

xlabel('Frequency (Hz)')
ylabel('Power')

xlim([10 20])


%% Using Matlab built-in function pwelch

srate = 1000; % in Hz
dt = 1/srate; % in s
f = 15; % in Hz
Tmax = 1;
t = dt:dt:Tmax;

X = sin(2*pi*f*t);

% subplot(3,1,1)
% plot(t,X)

window = length(X);
overlap = [];


% defining a frequency vector
freq = 0:0.01:80;

% freq = [0,1:5,6:0.5:14,14:0.01:20,25:5:80];
[Pxx,F]=pwelch(X,window,overlap,freq,srate);

% nfft defaults to the smaller x that satisfies: 2^x > window size

hold on
plot(F,Pxx,'y-','linew',3)

xlabel('Frequency (Hz)')
ylabel('Power (mV^2/Hz)')
hold off
% xlim([0 20])


%% Time-frequency uncertainty using pwelch


srate = 1000; % in Hz
dt = 1/srate; % in s
f = 5; % in Hz
Tmax = 10;
t = dt:dt:Tmax;

X = sin(2*pi*f*t);

% window = length(X)

window = 4*srate;
overlap = [];
% overlap  = window/2

nfft = 10000
freqresolution = srate/nfft

[Pxx,F]=pwelch(X,window,overlap,nfft,srate);

clf

plot(F,Pxx,'k-','linew',1)

xlabel('Frequency (Hz)')
ylabel('Power (mV^2/Hz)')

xlim([0 20])








