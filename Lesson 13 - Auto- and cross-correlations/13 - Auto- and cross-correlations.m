% Lesson 13

% Auto-correlation and cross-correlation

clear, clc, clf

X = ones(1,10);
Y = ones(1,3);


[CCG lags] = xcorr(X,Y)


% the third function argument (optional) defines the maximum analyzed lag
% [CCG lags] = xcorr(X,Y,10) 
[CCG lags] = xcorr(Y,X,10) 

% note that lag values are related to the 2nd time series


plot(lags,CCG,'ko-')

xlabel('Lag (ms)')
ylabel('CCG')

ylim([0 max(CCG)])

% autocorrelogram:

[ACG lags] = xcorr(X,X,10)

% ACG can be used with one single argument

[ACG lags] = xcorr(Y,10)

plot(lags,ACG,'ko-')

xlabel('Lag (ms)')
ylabel('ACG')

ylim([0 max(ACG)])


%% Applications...

% CCG can be used to:

% 1) assess synchrony
% 2) assess causality
% 3) assess oscillatory relationships

X = randn(1,3000);
% Y = randn(1,3000);

Y = circshift(X,124);

h1 = subplot(211);
plot(X)
xlabel('Time (ms)')

h2 = subplot(212);
plot(Y)
xlabel('Time (ms)')

linkaxes([h1 h2],'x')

% remember that the reference time series is the one passed as the second
% argument in 'xcorr'

clf

[CCG lags] = xcorr(Y,X);

plot(lags,CCG)

xlabel('Lag (ms)')
ylabel('CCG')

% ylim([-3000 3000])


%% Another example

clear, clc, clf

Stimulus = zeros(1,1000);

Istim = randi(975,[1,20]);

Stimulus(Istim) = 1;

Y = randi(40,[1,1000]);

Y(Istim+20) = 0*Y(Istim+20);
Y(Istim+21) = 0*Y(Istim+21);
Y(Istim+22) = 0*Y(Istim+22);
Y(Istim+23) = 0*Y(Istim+23);
Y(Istim+24) = 0*Y(Istim+24);

% Y(Istim+20) = 2*Y(Istim+20);
% Y(Istim+21) = 2*Y(Istim+21);
% Y(Istim+22) = 2*Y(Istim+22);
% Y(Istim+23) = 2*Y(Istim+23);
% Y(Istim+24) = 2*Y(Istim+24);

h1 = subplot(211);
plot(Y)
xlabel('Time (ms)')

h2 = subplot(212);
plot(Stimulus)
xlabel('Time (ms)')

linkaxes([h1 h2],'x')

%%
[CCG lags] = xcorr(Y,Stimulus,100);

clf
plot(lags,CCG)
xlabel('Lag (ms)')
ylabel('CCG')

%% Oscillatory example

clear

RandomSpikes = randi(1000,[1,100]);
SpikeTime(RandomSpikes) = 1;
SpikeTime(20:20:1000) = 1;

[ACG lags] = xcorr(SpikeTime,200);
ACG(lags==0) = 0;
plot(lags,ACG)

xlabel('Lag (ms)')
ylabel('ACG')

%% Oscillatory relationship between two signals
% i.e., cross-correlation

clear

srate = 1000;
dt = 1/srate;
t = dt:dt:3;

LFP1 = 1*sin(2*pi*4*t)+0*sin(2*pi*30*t)+0*randn(size(t));
LFP2 = 0*sin(2*pi*4*t+pi)+1*sin(2*pi*77*t)+0*randn(size(t));

subplot(211)
plot(t,LFP1)
hold on
plot(t,LFP2)
hold off
% xlim([1 1.5])

subplot(212)
[CCG lags] = xcorr(LFP2,LFP1)
plot(lags*dt,CCG)

xlabel('Time (s)')
ylabel('CCG')

% ylim([-1500 1500])

%% Application with real data
clear

load('LFP_HG_HFO.mat')

srate = 1000;
dt = 1/srate;
t = dt*(1:length(lfpHG));


% [ACG lags] = xcorr(lfpHG,1000);
[CCG lags] = xcorr(lfpHG,lfpHFO,1000);

subplot(211)
plot(t,lfpHG)
xlim([3 5])

subplot(212)
% plot(lags*dt,ACG)
plot(lags*dt,CCG)


%% CCG normalization

% To obtain the Pearson correlation coefficients, you can simply compute 
% the CCG between the z-scored time series

clear

srate = 1000;
dt = 1/srate;
t = dt:dt:3;

LFP1 = 1*sin(2*pi*4*t)+1*sin(2*pi*30*t)+0*randn(size(t));
LFP2 = 1*sin(2*pi*4*t+0*pi)+1*sin(2*pi*77*t)+0*randn(size(t));

normLFP1 = zscore(LFP1);
normLFP2 = zscore(LFP2);

clf

% [CCG lags] = xcorr(LFP2,LFP1)

[CCG lags] = xcorr(normLFP2,normLFP1)
CCG_pearson = CCG/(length(LFP1)-1);
plot(lags*dt,CCG_pearson)

xlabel('Time (s)')
ylabel('Pearson correlation')

% [CCG lags] = xcorr(normLFP2/sqrt(length(normLFP2)),...
%     normLFP1/sqrt(length(normLFP1)))
% 
% plot(lags*dt,CCG)
% xlabel('Time (s)')
% ylabel('Pearson correlation')

ylim([-1 1])

%%

[r p]=corr(LFP1(66:end)',LFP2(1:end-65)')

%% Matlab built-in CCG normalization

[CCGB lags] = xcorr(LFP2,LFP1,'coeff')

hold on

plot(lags*dt,CCGB)

hold off


%% Relationship between Power and ACG

LFP = sin(2*pi*8*t);

[ACG lags] = xcorr(LFP);

subplot(211)
plot(lags*dt,ACG)

%  The Fourier transform of a signal's ACG in a given frequency is 
%  equivalent to that frequency's Power

for f=1:20
   Power(f) = sum(ACG.*exp(-i*2*pi*f*(lags*dt))); 
end
    
subplot(212)
plot(1:20,real(Power),'bo-')
hold on
plot(1:20,imag(Power),'ro-')
hold off

























