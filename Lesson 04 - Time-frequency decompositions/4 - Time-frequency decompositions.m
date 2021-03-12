% Lesson 4

% Remember to look up the raw signal before doing any computations

clear 
clf
clc

load('LFP_HG_HFO.mat')

%%

srate = 1000; % in Hz (1/s)
dt = 1/srate; % in s
Tmax = length(lfpHG)*dt; % 300 s (5 min)
tvector = dt:dt:Tmax;
% tvector = (1:length(lfpHG))*dt;

plot(tvector,lfpHG)
hold on
plot(tvector,lfpHFO-1.2)
hold off
box off
xlim([0 1])


stepsize = 0.2
for nstep = 0:100

xlim([0 1]+nstep*stepsize)
   
ylim([-2 1])
pause(0.1)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Time-Frequency Power Decomposition (TFD)

clear, clf, clc

srate = 1000;
dt = 1/srate;
t = dt:dt:5;

LFP = 1*sin(2*pi*8*t);
LFP(1:2500) = 0;
LFP = LFP + 0.7*sin(2*pi*20*t);

subplot(311)
plot(t,LFP)
xlabel('time (s)')

win = 0.5*srate;
overlap = [];
nfft = 10000;

[Pxx F] = pwelch(LFP,win,overlap,nfft,srate);

subplot(3,1,[2 3])
plot(F,Pxx,'ko-')

xlabel('Frequency (Hz)')
ylabel('Power')
xlim([0 30])

%% TFD (this time for sure)

window_length = 0.5*srate;
step_size = 0.1*window_length; % (ie., 90% overlap)

Nwindows = (length(LFP)-window_length)/step_size+1;

clear T TFD

for nwin = 1:Nwindows
    winidx = (1:window_length) + (nwin-1)*step_size;
[Pxx F] = pwelch(LFP(winidx),window_length,[],nfft,srate);
T(nwin) = t(winidx(window_length/2));
TFD(nwin,:) = Pxx; 
    
end

% %

imagesc(T,F,TFD')
axis xy
ylim([0 30])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

% colorbar

%%

for nwin = 1:Nwindows
   plot(F,TFD(nwin,:))
   xlim([0 30])
    title(['Time = ' num2str(T(nwin)) ' s'])
   pause(1)
end
    

%% Using Matlab built-in function
 
window_length = 2*srate;
% step_size = 0.1*window_length; % (ie., 90% overlap)
overlap = 0.9*window_length;
nfft = 10000;

[S,F,T,TFD]=spectrogram(LFP,window_length,overlap,nfft,srate);

imagesc(T,F,TFD)

axis xy
ylim([0 30])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

colorbar


%% Normalizations

LFP = sin(2*pi*8*t);
LFP2 = sin(2*pi*8*t)+sin(2*pi*15*t);


[Pxx F] = pwelch(LFP,2*srate,[],nfft,srate);
[Pxx2 F] = pwelch(LFP2,2*srate,[],nfft,srate);

subplot(2,1,1)
plot(F,Pxx,'ko-')
hold on
plot(F,Pxx2,'ro-')
hold off
xlim([0 20])

ylabel('Power (mV^2/Hz)')
xlabel('Frequency (Hz)')

% Percentage power
NormPxx = Pxx/sum(Pxx);
NormPxx2 = Pxx2/sum(Pxx2);

subplot(2,1,2)
plot(F,NormPxx,'ko-')
hold on
plot(F,NormPxx2,'ro-')
hold off
xlim([0 20])

ylabel('Normalized Power (%)')
xlabel('Frequency (Hz)')


%% TFD normalizations

clear 

srate = 1000;
dt = 1/srate;
t=dt:dt:5;

LFP = sin(2*pi*8*t);

LFP(2500:end) =  LFP(2500:end)+0.1*sin(2*pi*50*t(2500:end));

window_length = 1*srate;
overlap = 0.9*window_length;
nfft = 10000;

[S,F,T,TFD]=spectrogram(LFP,window_length,overlap,nfft,srate);

subplot(211)

imagesc(T-2.5,F,TFD)
hold on

plot([0 0],[0 100],'w--','linew',2)
hold off
axis xy
ylim([0 100])
xlabel('Time (s)')
ylabel('Frequency (Hz)')


% Baseline normalization

I = find(T<2.5)

PSDbaseline = mean(TFD(:,I),2);

NormTFD = TFD./repmat(PSDbaseline,1,length(T));

subplot(212)

imagesc(T-2.5,F,NormTFD)
hold on

plot([0 0],[0 100],'w--','linew',2)
hold off
axis xy
ylim([0 100])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

















