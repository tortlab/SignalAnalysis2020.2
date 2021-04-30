% Lesson 14

% Cross-frequency coupling

% Amplitude
% Phase
% Frequency

% Phase-Phase (n:m phase-locking)
% Phase-Amplitude
% Amplitude-Amplitude
% Phase-Frequency
% Amplitude-Frequency
% Phase-Frequency
% Frequency-Frequency

%% n:m phase-locking

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:5;

w = 2*pi*8 + 5*randn(size(t));
theta = sin(cumsum(w*dt)+pi/2);

w = 2*pi*12 + 50*randn(size(t));
gamma = sin(cumsum(w*dt));

thetaphase = angle(hilbert(theta));
gammaphase = angle(hilbert(gamma));

n = 1
m = 5

Mthetaphase = angle(exp(1i*m*thetaphase));
Ngammaphase = angle(exp(1i*n*gammaphase));


DeltaPhase = angle(exp(1i*(thetaphase-gammaphase)));
NMDeltaPhase = angle(exp(1i*(m*thetaphase-n*gammaphase)));



subplot(311)
plot(t,theta)
hold on
plot(t,gamma)
hold off
xlim([0 .5])

subplot(312)
plot(t,Mthetaphase,'b.')
hold on
plot(t,Ngammaphase,'r.')
hold off
xlim([0 .5])
ylim([-pi pi])

subplot(313)
plot(t,NMDeltaPhase,'k.')
ylim([-pi pi])
xlim([0 .5])

% PLV = abs(mean(exp(i*DeltaPhase)));
% rose(DeltaPhase,18)
% title(['Phase-Locking Value = ' num2str(PLV)])


NMPLV = abs(mean(exp(i*NMDeltaPhase)));
rose(NMDeltaPhase,18)
title(['N:M Phase-Locking Value = ' num2str(NMPLV)])

%% Finding best "m"

n=1

for m=1:10
    
NMDeltaPhase = angle(exp(1i*(m*thetaphase-n*gammaphase)));    
NMPLV = abs(mean(exp(i*NMDeltaPhase)));

nmPLV(m) = NMPLV;
end

clf
plot(1:10,nmPLV,'ko-')
set(gca,'xtick',1:10)
xlabel('m')
ylabel('N:M PLV')
ylim([0 1])

%% n:m phase locking comodulogram

clear nmPLV

for n=1:10
for m=1:10
    
NMDeltaPhase = angle(exp(1i*(m*thetaphase-n*gammaphase)));    
NMPLV = abs(mean(exp(i*NMDeltaPhase)));
nmPLV(n,m) = NMPLV;   

end
end

clf
imagesc(nmPLV')


axis xy
xlabel('n')
ylabel('m')
axis square

%% Phase-amplitude coupling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clf, clc

load('LFP_HG_HFO.mat')

srate = 1000;
dt = 1/srate;
t = dt*(1:length(lfpHG));

LFP = lfpHG;

phase_freq = eegfilt(LFP,srate,6,10);
amp_freq = eegfilt(LFP,srate,60,90);

phase = angle(hilbert(phase_freq));
ampenv = abs(hilbert(amp_freq));

figure(1)

h1 = subplot(211)
plot(t,LFP,'k-')
hold on
plot(t,phase_freq-1,'b-','linew',2)
plot(t,5*amp_freq-2,'r-','linew',0.5)
plot(t,5*ampenv-2,'r-','linew',2)

hold off

h2 = subplot(212)
plot(t,phase,'b.')

linkaxes([h1 h2],'x')

xlim([3 4])

% % identifying phase bins

% I = find(phase > deg2rad(20) & phase < deg2rad(40));

j=17
I = find(phase > deg2rad(-180+j*20) & phase < deg2rad(-160+j*20));


hold on
plot(t(I),phase(I),'y.')

subplot(211)
hold on
plot(t(I),5*ampenv(I)-2,'y.')
hold off

MeanAmp = mean(ampenv(I))
title(['Mean \gamma Amplitude = ' num2str(MeanAmp)])


% % looping over phase bins
clear MeanAmp

count = 0;

for phasebin = -180:20:160
count = count+1;
I = find(phase > deg2rad(phasebin) & phase < deg2rad(phasebin+20));
MeanAmp(count) = mean(ampenv(I));
end

% %

figure(2)
% subplot(121)
% bar(10:20:710,[MeanAmp MeanAmp])
% xlabel('\Theta Phase (Deg)')
% ylabel('Mean \gamma Amplitude')
% 
% set(gca,'xtick',0:90:720)


% normalizing mean gamma amplitude in order to obtain a probability
% distribution, that is, a distribution that satisfies pJ>=0 and
% sum(pJ)=0

p = MeanAmp/sum(MeanAmp);

subplot(111)
bar(10:20:710,[p p])
xlabel('\Theta Phase (Deg)')
ylabel('Normalized mean \gamma Amplitude')

set(gca,'xtick',0:90:720)

% computing entropy:

H = -sum(p.*log(p));
Hmax = log(length(p));
MI = (Hmax-H)/Hmax;

title(['Modulation Index = ' num2str(MI)])

%% Comodulograma



clear, clf, clc

%%
load('LFP_HG_HFO.mat')

srate = 1000;
dt = 1/srate;
t = dt*(1:length(lfpHG));

LFP = lfpHFO;

clear Comodulogram

phase_freq_vector = 0:2:20;
fp_bandwidth = 4;

amp_freq_vector = 20:5:200;
fa_bandwidth = 20;

count_phase = 0;
tic
for fp = phase_freq_vector
count_phase = count_phase+1

phase_freq = eegfilt(LFP,srate,fp,fp+fp_bandwidth);
phase = angle(hilbert(phase_freq));

count_amp = 0;
for fa = amp_freq_vector
count_amp = count_amp+1;

amp_freq = eegfilt(LFP,srate,fa,fa+fa_bandwidth);
ampenv = abs(hilbert(amp_freq));

% MI calculation
count = 0;
for phasebin = -180:20:160
count = count+1;
I = find(phase > deg2rad(phasebin) & phase < deg2rad(phasebin+20));
MeanAmp(count) = mean(ampenv(I));
end

p = MeanAmp/sum(MeanAmp);

H = -sum(p.*log(p));
Hmax = log(length(p));
MI = (Hmax-H)/Hmax;

Comodulogram(count_phase,count_amp)=MI;

end
end
toc

%%

% imagesc(phase_freq_vector+fp_bandwidth/2,...
%     amp_freq_vector+fa_bandwidth/2,...
%     Comodulogram')
% axis xy
% xlabel('Phase Frequency (Hz)')
% ylabel('Amplitude Frequency (Hz)')

subplot(122)
contourf(phase_freq_vector+fp_bandwidth/2,...
    amp_freq_vector+fa_bandwidth/2,...
    Comodulogram',30,'linestyle','none')
axis xy
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')

%% Alternative method, with less filtering

tic
PhaseAll = zeros(length(phase_freq_vector),length(LFP));
count_phase = 0;
for fp = phase_freq_vector
count_phase = count_phase+1

phase_freq = eegfilt(LFP,srate,fp,fp+fp_bandwidth);
phase = angle(hilbert(phase_freq));
PhaseAll(count_phase,:)=phase;
end

AmpAll = zeros(length(amp_freq_vector),length(LFP));
count_amp = 0;
for fa = amp_freq_vector
count_amp = count_amp+1;

amp_freq = eegfilt(LFP,srate,fa,fa+fa_bandwidth);
ampenv = abs(hilbert(amp_freq));
AmpAll(count_amp,:)=ampenv;
end

count_phase = 0;
for fp = phase_freq_vector
count_phase = count_phase+1

count_amp = 0;
for fa = amp_freq_vector
count_amp = count_amp+1;

% MI calculation
count = 0;
for phasebin = -180:20:160
count = count+1;
I = find(PhaseAll(count_phase,:) > deg2rad(phasebin)...
    & PhaseAll(count_phase,:)  < deg2rad(phasebin+20));
MeanAmp(count) = mean(AmpAll(count_amp,I));
end

p = MeanAmp/sum(MeanAmp);

H = -sum(p.*log(p));
Hmax = log(length(p));
MI = (Hmax-H)/Hmax;

Comodulogram(count_phase,count_amp)=MI;

end
end
toc


%% Althernative method 2

tic
PhaseAll = zeros(length(phase_freq_vector),length(LFP));
count_phase = 0;
for fp = phase_freq_vector
count_phase = count_phase+1

phase_freq = eegfilt(LFP,srate,fp,fp+fp_bandwidth);
phase = angle(hilbert(phase_freq));
PhaseAll(count_phase,:)=phase;
end

count_amp = 0;
for fa = amp_freq_vector
count_amp = count_amp+1;

amp_freq = eegfilt(LFP,srate,fa,fa+fa_bandwidth);
ampenv = abs(hilbert(amp_freq));

count_phase = 0;
for fp = phase_freq_vector
count_phase = count_phase+1;

% MI calculation
count = 0;
for phasebin = -180:20:160
count = count+1;
I = find(PhaseAll(count_phase,:) > deg2rad(phasebin)...
    & PhaseAll(count_phase,:)  < deg2rad(phasebin+20));
MeanAmp(count) = mean(ampenv(I));
end

p = MeanAmp/sum(MeanAmp);

H = -sum(p.*log(p));
Hmax = log(length(p));
MI = (Hmax-H)/Hmax;

Comodulogram(count_phase,count_amp)=MI;
end
end
toc

%% Phase-Energy Plot

load('LFP_HG_HFO.mat')

srate = 1000;
dt = 1/srate;
t = dt*(1:length(lfpHG));

LFP = lfpHFO;

theta = eegfilt(LFP,srate,5,10);

amp_freq_vector = 20:5:200;
fa_bandwidth = 20;

AmpAll = zeros(length(amp_freq_vector),length(LFP));
count_amp = 0;
for fa = amp_freq_vector
count_amp = count_amp+1;

amp_freq = eegfilt(LFP,srate,fa,fa+fa_bandwidth);
ampenv = abs(hilbert(amp_freq));
AmpAll(count_amp,:)=ampenv;
end

% %
clf

plot(t,theta)
xlim([50 55])

[pks I] = findpeaks(theta);

hold on
plot(t(I),pks,'*')
plot(t(I),theta(I),'bo')
hold off


windowlength = 1*srate;

EnergyPhase = zeros(length(amp_freq_vector),windowlength+1);
meanLFP = zeros(1,windowlength+1);


count = 0
for j=1:length(I)

if I(j) > windowlength/2 & I(j) < length(LFP)-windowlength/2
    count=count+1
    WindowIndex = I(j)-windowlength/2:I(j)+windowlength/2;
    
  EnergyPhase = EnergyPhase  + AmpAll(:,WindowIndex);
  meanLFP =  meanLFP + LFP(WindowIndex);
end

end

EnergyPhase = EnergyPhase/count;
meanLFP = meanLFP/count;

% %

clear EnergyPhaseNorm
for j = 1:length(amp_freq_vector)
  EnergyPhaseNorm(j,:) = EnergyPhase(j,:)/mean(EnergyPhase(j,:));  
end

subplot(5,1,[1 4])
contourf(-0.5:dt:0.5,...
    amp_freq_vector+fa_bandwidth/2,...
    EnergyPhaseNorm,30,'linestyle','none')
axis xy
ylabel('Frequency (Hz)')

subplot(515)
plot(-0.5:dt:0.5,meanLFP,'k-','linew',2)
xlabel('time (s)')
