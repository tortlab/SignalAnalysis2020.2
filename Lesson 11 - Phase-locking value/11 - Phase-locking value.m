% Lesson 11

% Obtaining the instantaneous phase of the signal

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:2;

A = 2+sin(2*pi*1*t);
LFP = A.*sin(2*pi*7*t);

AnalyticLFP = hilbert(LFP);
HT = imag(AnalyticLFP);
EnvAmp = abs(AnalyticLFP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Phase = angle(AnalyticLFP); % radians from -pi to pi
% Phase = Phase+pi;
% Phase(Phase<0) = Phase(Phase<0)+2*pi;

Phase = rad2deg(Phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1:2:2000

h1 = subplot(311)
plot(t,LFP,'b-')
hold on
plot(t,HT,'r-')
plot(t,EnvAmp,'k','linew',2)

plot(t(j),LFP(j),'bo','markerf','b')
plot(t(j),HT(j),'ro','markerf','r')
plot(t(j),EnvAmp(j),'ko','markerf','k')
hold off
xlabel('Time (s)')
ylabel('mV')

h2 = subplot(312)
plot(t,Phase,'g.')
hold on
plot(t(j),Phase(j),'go','markerf','g')

hold off
% ylim([0 2*pi])
% ylim([0 360])
% set(gca,'ytick',0:90:360)
set(gca,'ytick',-180:90:180)
title(['Inst Phase = ' num2str(Phase(j)) '^o'])

xlabel('Time (s)')
ylabel('Phase (^o)')


subplot(313)
plot([0 LFP(j)],[0 HT(j)],'k-')
hold on
plot(LFP(j),HT(j),'ko','markerf','k')
plot(LFP(j),0,'bo','markerf','b')
plot(0,HT(j),'ro','markerf','r')
hold off

xlim([-3 3])
ylim([-3 3])
axis square

pause(0.001)
end

linkaxes([h1 h2],'x')

%% Instantaneous frequency

t = dt:dt:4;

A = 2+sin(2*pi*1*t);
LFP = A.*sin(2*pi*7*t);

% Frequency modulation

% F = 2*pi*10*t;
% freq = 10+4*sin(2*pi*1*t);
% freq = 10 + randn(size(t));
% LFP = A.*sin(2*pi*freq.*t);

% programming frequency modulation properly
freq = 10+4*sin(2*pi*1*t);
LFP = A.*sin(cumsum(2*pi*freq.*dt));

AnalyticLFP = hilbert(LFP);
HT = imag(AnalyticLFP);
EnvAmp = abs(AnalyticLFP);

Phase = angle(hilbert(LFP));
% % % Freq = diff(Phase)/(2*pi*dt); % wrong equation
PhaseUnwrapped = unwrap(Phase);
Freq = diff(PhaseUnwrapped)/(2*pi*dt); 

% alternative notation
Freq2 = angle(exp(1i*diff(Phase)))/(2*pi*dt);


h1 = subplot(311);
plot(t,LFP,'b-')
hold on
plot(t,HT,'r-')
plot(t,EnvAmp,'k','linew',2)
hold off
xlabel('Time (s)')
ylabel('mV')

h2 = subplot(312);
plot(t,Phase,'g.')
% hold on
% plot(t,PhaseUnwrapped,'g-')
% hold off
xlabel('Time (s)')
ylabel('Phase (rad)')

h3 = subplot(313);
plot(t(1:end-1)+dt/2,Freq)
hold on
plot(t(1:end-1)+dt/2,Freq2)
hold off
% ylim([0 15])


linkaxes([h1 h2 h3],'x')
xlabel('Time (s)')
ylabel('Inst Freq (Hz)')

%% Phase-locking metrics

clear, clc, clf

srate = 1000;
dt = 1/srate;
t = dt:dt:4;

noise1 = 3*randn(size(t));
noise2 = 3*randn(size(t));

LFP1 = sin(2*pi*8*t) + noise1;
LFP2 = sin(2*pi*8*t+pi/2) + noise2;

LFP1filtered = eegfilt(LFP1,srate,6,10);
LFP2filtered = eegfilt(LFP2,srate,6,10);

% LFP1filtered = eegfilt(LFP1,srate,50,60);
% LFP2filtered = eegfilt(LFP2,srate,50,60);


Phase1 = angle(hilbert(LFP1filtered));
Phase2 = angle(hilbert(LFP2filtered));

% DeltaPhase = Phase2-Phase1; % > can be problematic
DeltaPhase = angle(exp(1i*(Phase2-Phase1))); % recommended method

figure(1)

h1 = subplot(311)
plot(t,LFP1,'b')
hold on
plot(t,LFP2,'r')
plot(t,LFP1filtered-3,'b')
plot(t,LFP2filtered-3,'r')
hold off

xlim([2 3])

h2 = subplot(312)
plot(t,Phase1,'b.')
hold on
plot(t,Phase2,'r.')
hold off
ylabel('\Phi (rad)')
set(gca,'ytick',-pi:pi/2:pi,...
    'yticklabel',{'-pi','-pi/2','0','pi/2','pi'})
xlim([2 3])


h3 = subplot(313)
plot(t,rad2deg(DeltaPhase),'k.')
ylabel('\Delta\Phi (^o)')
xlabel('Time (s)')
xlim([2 3])
ylim([-180 180])

linkaxes([h1 h2 h3],'x')

% % phase difference histograms

figure(2)

subplot(221)

phasebins = -170:20:170

[counts phasebins]=hist(DeltaPhase,deg2rad(phasebins));

bar(rad2deg(phasebins),counts,'k')
xlabel('\Delta\Phi (^o)')
ylabel('counts')
set(gca,'xtick',-180:90:180)

subplot(222)
rose(DeltaPhase,18) % number of phase bins
title('\Delta\Phi (^o)')

subplot(223)
% z1 = 1 + 1*1i;
% z2 = -1 + 1*1i;
% compass([z1 z2]);

compass(exp(1i*DeltaPhase(1:20:end)))


subplot(224)

polar(angle(exp(1i*DeltaPhase(1:20:end))),...
    abs(exp(1i*DeltaPhase(1:20:end))),'ko')

MeanVector = mean(exp(1i*DeltaPhase));

hold on
compass(MeanVector,'r')
hold off


% computing the Phase-Locking Value
PLV = abs(mean(exp(1i*DeltaPhase)));

title(['PLV = ' num2str(PLV)])











