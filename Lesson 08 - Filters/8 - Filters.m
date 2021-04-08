% Lesson 08

% Filters using eegfilt

clear
clc
clf

srate = 1000;
dt = 1/srate;
t = dt:dt:2;

LFP = sin(2*pi*10*t)+sin(2*pi*30*t)+sin(2*pi*50*t);

% eegfilt calls the function filtfilt, which applies a filter
% two times (the 2nd is on reverse, to correct for phase distortions)

% default filter order is given by the formula:
% default = 3*fix(srate/locutoff)

% Low-pass filter
% filtered = eegfilt(LFP,srate,0,16);

% High-pass filter
% filtered = eegfilt(LFP,srate,35,0);

% Band-pass filter
% filtered = eegfilt(LFP,srate,45,55);

% changing filter order
% order = 66
% filtered = eegfilt(LFP,srate,45,55,0,order);

% Notch filter
filtered = eegfilt(LFP,srate,20,40,0,[],1);


[Pxx F] = pwelch(LFP,length(LFP),[],2^16,srate);
[PxxFilt F] = pwelch(filtered,length(LFP),[],2^16,srate);

figure(2)

subplot(311)
plot(t,LFP)
hold on
plot(t,filtered,'r-')
hold off

subplot(312)
plot(F,Pxx)
hold on
plot(F,PxxFilt,'r-')
hold off
xlim([0 70])

subplot(313)
plot(F,PxxFilt,'r-')
xlim([0 70])


%% No filter is perfect

figure(1)


LowFreqCutoff = 55;
HighFreqCutoff = 65;

whitenoise = randn(1,1000000);
 
% filtered = eegfilt(whitenoise,srate,...
%     LowFreqCutoff,HighFreqCutoff);
% 
% order = 100
% filtered = eegfilt(whitenoise,srate,...
%     LowFreqCutoff,HighFreqCutoff,0,order);

filtered = eegfilt(whitenoise,srate,...
    LowFreqCutoff,HighFreqCutoff,0,[],1);

[PxxW F] = pwelch(whitenoise,srate,[],2^16,srate);
[PxxF F] = pwelch(filtered,srate,[],2^16,srate);

subplot(211)
plot(F,PxxW,'k-')
hold on
plot(F,PxxF)
plot([HighFreqCutoff HighFreqCutoff],[0 max(PxxW)*1.2],'k-')
plot([LowFreqCutoff LowFreqCutoff],[0 max(PxxW)*1.2],'k-')

hold off
xlim([0 100])

%% understanding eegfilt's 5th parameter

% epochframes = frames per epoch 
%(filter each epoch separately {def/0: data is 1 epoch}

t = dt:dt:1;

LFP1 = sin(2*pi*10*t);
LFP2 = sin(2*pi*10*t+pi);
LFP3 = sin(2*pi*10*t+pi/2);
LFP4 = sin(2*pi*10*t+pi/3);

LFP = [LFP1 LFP2 LFP3 LFP4];

filteredTrials = eegfilt(LFP,srate,0,15,1000);
filteredSingle = eegfilt(LFP,srate,0,15);

clf

plot(LFP,'k-')
hold on
plot(filteredSingle)
plot(filteredTrials,'r-')
hold off

%% Order effects again (temporal leakage)

srate = 1000
dt = 1/srate
t = dt:dt:4;

lfp = sin(2*pi*10*t);
lfp(1:2000)=0;
lfp(3000:end)=0;

order = 500
filtered = eegfilt(lfp,srate,0,15,0,order);

plot(t,lfp)
hold on
plot(t,filtered,'r-','linew',2)
hold off


xlim([1.5 3.5])

%% No filter is perfect (again)


srate = 1000;
dt = 1/srate;
t = dt:dt:4;

lfp = sin(2*pi*10*t);
lfp = lfp + 0.0*randn(size(t));

filtered = eegfilt(lfp,srate,15,0);

[Pxx F] = pwelch(lfp,4*srate,[],2^16,srate);
[PxxF F] = pwelch(filtered,4*srate,[],2^16,srate);


subplot(221)
plot(t,lfp)
hold on
plot(t,filtered)
hold off
xlim([1 1.5])

subplot(222)
plot(F,Pxx)
hold on
plot(F,PxxF)
hold off
xlim([0 20])

subplot(223)
plot(t,filtered,'r-')
xlim([1 1.5])

subplot(224)
plot(F,PxxF,'r-')
xlim([0 20])

%% eegfilt accepts matrices as inputs

% in this case it applies filters line-by-line

% Channel vs timepoints

LFPAll(1,:) = sin(2*pi*5*t);
LFPAll(2,:) = sin(2*pi*15*t);
LFPAll(3,:) = sin(2*pi*30*t);

filtered = eegfilt(LFPAll,srate,10,25,0,[],1);


h1 = subplot(311);
plot(t,LFPAll(1,:))
hold on
plot(t,filtered(1,:))
hold off
ylim([-3 3])

h2 = subplot(312);
plot(t,LFPAll(2,:))
hold on
plot(t,filtered(2,:))
hold off
ylim([-3 3])

h3 = subplot(313);
plot(t,LFPAll(3,:))
hold on
plot(t,filtered(3,:))
hold off
ylim([-3 3])

linkaxes([h1 h2 h3])

% linkaxes([h1 h2 h3],'x')