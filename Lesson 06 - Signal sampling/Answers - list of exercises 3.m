% Answers to list of exercises 3 (first assessment)

load('LFPprobe.mat')

dt = 1/srate;
t = (1:length(LFPprobe))*dt;

%% Question 1

clf

% Plot all channels simultaneously (without overlapping)

size(LFPprobe)

for ch = 1:16
plot(t,LFPprobe(ch,:)-ch*1000,'k-')
hold on
end
hold off

set(gca,'ytick',-16000:1000:-1000)
set(gca,'TickDir','out')
set(gca,'yticklabel',1600:-100:100)

box off

xlim([96 98])
xlabel('Time (s)')
ylabel('Depth (\mum)')

%% Alternative solution

shift = [-1000:-1000:-16000]';
LFPrescaled = LFPprobe + shift;
% for older Matlab versions:
% LFPrescaled = LFPprobe + repmat(shift,1,length(LFPprobe));
plot(t,LFPrescaled,'k-')

set(gca,'ytick',-16000:1000:-1000)
set(gca,'TickDir','out')
set(gca,'yticklabel',1600:-100:100)

box off

xlim([96 98])
xlabel('Time (s)')
ylabel('Depth (\mum)')



%% Question 2
% Vary the period shown with 2-s sliding window (LFPtv).

for nstep = 0:500
   xlim([0 2] + nstep*0.1)
   pause(0.0005)
end


%% Question 3

% Compute the PSD for all channels, using both absolute
% and normalized power

clear PxxAll NormPxxAll

PxxAll = zeros(16,2001);
NormPxxAll = zeros(16,2001);

tic
for ch = 1:16
   [Pxx F] = pwelch(LFPprobe(ch,:),4*srate,[],0:0.1:200,srate); 
   PxxAll(ch,:) = Pxx;
   NormPxxAll(ch,:) = Pxx/sum(Pxx);   
end
toc

%% alternative method without loop:

% When X is a matrix, the PSD is
% computed independently for each column and 
% stored in the corresponding column of Pxx.

tic
[PxxAll2 F] = pwelch(LFPprobe',4*srate,[],0:0.1:200,srate);
PxxAll2 = PxxAll2';
NormPxxAll2 = PxxAll2./sum(PxxAll2,2);
NormPxxAll2b = PxxAll2./repmat(sum(PxxAll2,2),1,2001);
toc

%% Question 4


% Plot four panels on different subplots, 
% where each panel shows the PSDs of all channels 
% in the following variations: 
% [absolute x normalized] X [linear x logarithmic scale]

for ch = 1:16
   Label{ch} = ['Ch ' int2str(ch)]; 
end


subplot(221)
plot(F,PxxAll)
xlabel('Freq (Hz)')
xlim([0 20])
ylabel('Power')
box off

% legend(Label)
% legend('boxoff')
% legend('location','northeastoutside')

subplot(222)
plot(F,NormPxxAll)
xlabel('Freq (Hz)')
xlim([0 20])
ylabel('Power')
box off

subplot(223)
plot(F,10*log10(PxxAll))
xlabel('Freq (Hz)')
xlim([0 200])
ylabel('Power (dB)')
box off


subplot(2,2,4)
I = find(F < 59 | F > 61);
plot(F(I),10*log10(NormPxxAll(:,I)))
xlabel('Freq (Hz)')
xlim([0 200])
ylabel('Power (dB)')
box off

%% Question 5

% Compute the average power (absolute and normalized) 
% on the theta frequency range (5-10 Hz) for each channel

I = find(F > 5 & F <10);

ThetaPower = mean(PxxAll(:,I),2);
NormThetaPower = mean(NormPxxAll(:,I),2);

%% Question 6

% Create a line plot of the values obtained in 5), 
% indicating the spatial location (‘height’) of the electrode 
% (from 100 to 1600?m) on the Y axis and average power on the X axis. 
% Do one figure for absolute values and another for normalized ones.

subplot(121)

plot(ThetaPower,-100:-100:-1600,'ko-',...
    'markersize',12,'markerfacecolor','y',...
    'linewidth',2)

set(gca,'ytick',-1600:100:-100,...
    'yticklabel',1600:-100:100)

xlabel('Theta Power')
ylabel('Depth (\mum)')

ylim([-1700 0])

subplot(122)

plot(NormThetaPower,-100:-100:-1600,'ko-',...
    'markersize',12,'markerfacecolor','y',...
    'linewidth',2)

set(gca,'ytick',-1600:100:-100,...
    'yticklabel',1600:-100:100)

xlabel('Norm Theta Power')
ylabel('Depth (\mum)')

ylim([-1700 0])

%% Question 7

% Repeat 5) and 6) for gamma (30 to 100 Hz) 
% and HFO (120 to 160 Hz) frequency bands

I = find(F > 60 & F <100);
GammaPower = mean(PxxAll(:,I),2);
NormGammaPower = mean(NormPxxAll(:,I),2);

% I = find(F > 120 & F <160);
% HFOPower = mean(PxxAll(:,I),2);
% NormHFOPower = mean(NormPxxAll(:,I),2);

subplot(121)

plot(GammaPower,-100:-100:-1600,'ko-',...
    'markersize',12,'markerfacecolor','y',...
    'linewidth',2)

set(gca,'ytick',-1600:100:-100,...
    'yticklabel',1600:-100:100)

xlabel('Gamma Power')
ylabel('Depth (\mum)')

ylim([-1700 0])

subplot(122)

plot(NormGammaPower,-100:-100:-1600,'ko-',...
    'markersize',12,'markerfacecolor','y',...
    'linewidth',2)

set(gca,'ytick',-1600:100:-100,...
    'yticklabel',1600:-100:100)

xlabel('Norm Gamma Power')
ylabel('Depth (\mum)')

ylim([-1700 0])

%% Question 8

% Compute and plot the TFD for channels 1 and 16 in different subplots

[S,F,T,TFDch1]=spectrogram(LFPprobe(1,:),4*srate,0.9*4*srate,2^15,srate);
[S,F,T,TFDch16]=spectrogram(LFPprobe(16,:),4*srate,0.9*4*srate,2^15,srate);

h1 = subplot(211);
imagesc(T,F,TFDch1)
axis xy
ylim([0 30])
xlabel('Time (s)')
ylabel('Freq (Hz)')
h = colorbar
ylabel(h,'Power','fontsize',12)

h2 = subplot(212);
imagesc(T,F,TFDch16)
axis xy
ylim([0 30])
xlabel('Time (s)')
ylabel('Freq (Hz)')
h = colorbar
ylabel(h,'Power','fontsize',12)

linkaxes([h1 h2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question 9

% Use the result from 8) to compute and plot the time series of 
% average power in theta, gamma and HFO on both channels
% (show each frequency band in a different subplot)

Itheta = find(F>5 & F<10);
Igamma = find(F>30 & F<100);
IHFO = find(F>120 & F<160);

Theta1 = mean(TFDch1(Itheta,:));
Theta16 = mean(TFDch16(Itheta,:));
Gamma1 = mean(TFDch1(Igamma,:));
Gamma16 = mean(TFDch16(Igamma,:));
HFO1 = mean(TFDch1(IHFO,:));
HFO16 = mean(TFDch16(IHFO,:));

subplot(311)
plot(T,Theta1)
hold on
plot(T,Theta16)
hold off
xlabel('Time')
ylabel('\theta power')

subplot(312)
plot(T,Gamma1)
hold on
plot(T,Gamma16)
hold off
xlabel('Time')
ylabel('\gamma power')

subplot(313)
plot(T,HFO1)
hold on
plot(T,HFO16)
hold off
xlabel('Time')
ylabel('HFO power')


%% alternative with normalized values

subplot(311)
plot(T,Theta1/mean(Theta1))
hold on
plot(T,Theta16/mean(Theta16))
hold off
xlabel('Time')
ylabel('\theta power')

subplot(312)
plot(T,Gamma1/mean(Gamma1))
hold on
plot(T,Gamma16/mean(Gamma16))
hold off
xlabel('Time')
ylabel('\gamma power')

subplot(313)
plot(T,HFO1/mean(HFO1))
hold on
plot(T,HFO16/mean(HFO16))
hold off
xlabel('Time')
ylabel('HFO power')

%% z-score normalization


subplot(311)
plot(T,zscore(Theta1))
hold on
plot(T,zscore(Theta16))
hold off
xlabel('Time')
ylabel('norm \theta power')

subplot(312)
plot(T,zscore(Gamma1))
hold on
plot(T,zscore(Gamma16))
hold off
xlabel('Time')
ylabel('norm \gamma power')

subplot(313)
plot(T,zscore(HFO1))
hold on
plot(T,zscore(HFO16))
hold off
xlabel('Time')
ylabel('norm HFO power')

%% Question 10

% Plot, in different subplots, scatterplots of 
% every possible combination of the time series in 9) 
% (e.x., theta ch 1 vs theta ch 16, theta ch 1 vs gamma ch 16, etc), 
% using their linear correlation value as the title of each subplot

% showing only first four combinations as proof of concept

subplot(221)
plot(Theta1,Theta16,'ko')
[r p] = corr(Theta1',Theta16');
title(['r = ' num2str(r)])
xlabel('Theta Ch1')
ylabel('Theta Ch16')
axis square
axis tight

subplot(222)
plot(Gamma1,Gamma16,'ko')
[r p] = corr(Gamma1',Gamma16');
title(['r = ' num2str(r)])
xlabel('Gamma Ch1')
ylabel('Gamma Ch16')
axis square
axis tight

subplot(223)
plot(HFO1,HFO16,'ko')
[r p] = corr(HFO1',HFO16');
title(['r = ' num2str(r)])
xlabel('HFO Ch1')
ylabel('HFO Ch16')
axis square
axis tight

subplot(224)
plot(Theta1,Gamma16,'ko')
[r p] = corr(Theta1',Gamma16');
title(['r = ' num2str(r)])
xlabel('Theta Ch1')
ylabel('Gamma Ch16')
axis square
axis tight

%% Question 11

clf

% Compute and plot the coherence spectrum of 
% channels 2 and 16 with channel 1, separately

clear Label CxyAll
for ch=2:16
[Cxy F] = mscohere(LFPprobe(1,:),LFPprobe(ch,:),4*srate,[],2^15,srate);

CxyAll(ch-1,:)=Cxy;

plot(F,Cxy)
hold on

Label{ch-1} = ['Ch1-Ch' num2str(ch)];
end
hold off
% %

xlim([0 30])

xlabel('Freq (Hz)')
ylabel('Coherence')
legend(Label,'location','northeastoutside')

set(gcf,'color','white')

%% Question 12

% Compute the average coherence on theta and gamma
% for the same channel pairs in 11)

MeanThetaCoh = mean(CxyAll(:,F>5 & F<10),2);
MeanGammaCoh = mean(CxyAll(:,F>30 & F<100),2);

%% Question 13

% Create a line plot of the values obtained in 12), 
% indicating the spatial location (‘height’) of the electrode 
% (from 100 to 1600?m) on the Y axis 
% and average coherence values on the X axis.

subplot(121)

plot(MeanThetaCoh,-200:-100:-1600,'ko-',...
    'markersize',12,'markerfacecolor','y',...
    'linewidth',2)

set(gca,'ytick',-1600:100:-200,...
    'yticklabel',1600:-100:200)

xlabel('Theta Coherence with Ch 1')
ylabel('Depth (\mum)')

ylim([-1700 -100])

subplot(122)

plot(MeanGammaCoh,-200:-100:-1600,'ko-',...
    'markersize',12,'markerfacecolor','y',...
    'linewidth',2)

set(gca,'ytick',-1600:100:-200,...
    'yticklabel',1600:-100:200)

xlabel('Gamma Coherence with Ch 1')
ylabel('Depth (\mum)')

ylim([-1700 -100])

%% Question 14

% Compute and plot the coherogram between channels 1 and 16

win = 20*srate; % window for Cxy computation
step = 0.1*win;
Nwin = (length(t)-win)/step+1;

cohwin = 2*srate;
nfft = 2^16;
overlap = 0;

% X = LFPprobe(1,:);
% Y = LFPprobe(16,:);

for nwin = 1:Nwin
    nwin
 winidx = (1:win) + (nwin-1)*step; 
%  [Cxy F] = mscohere(X(winidx),Y(winidx),cohwin,overlap,nfft,srate);
  
 [Cxy F] = mscohere(LFPprobe(1,winidx),...
     LFPprobe(16,winidx),cohwin,overlap,nfft,srate);
 Coherogram(nwin,:) = Cxy;
  T(nwin) = t(winidx(win/2));    
end
    
% %
clf
imagesc(T,F,Coherogram')
axis xy
ylim([0 20])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

h = colorbar
h.Label.String = 'Coherence'

%% Question 15

% Compute and plot (using imagesc) the average coherence in
% a given frequency band (e.x., theta or gamma) 
% for every possible pair of channels.

clear MeanThetaCohAll MeanGammaCohAll

% for ch1 = 1:16
%    for ch2 = 1:16
for ch1 = 1:15
   for ch2 = ch1+1:16
              
       [ch1 ch2]
      [Cxy F] = mscohere(LFPprobe(ch1,:),LFPprobe(ch2,:),...
          10*srate,[],2^15,srate); 
       
 MeanThetaCohAll(ch1,ch2-1) = mean(Cxy(F>5 & F<10));
MeanGammaCohAll(ch1,ch2-1) = mean(Cxy(F>30 & F<100));   
   end
end

% %

% ThetaCoh2 = triu(MeanThetaCohAll) - eye(16);
% 
% imagesc(1:16,1:16,MeanThetaCohAll)
% 
% imagesc(1:16,1:16,ThetaCoh2)

imagesc(2:16,1:15,MeanThetaCohAll)
set(gca,'ytick',1:15,'xtick',2:16)

xlabel('Ch #')
ylabel('Ch #')

axis square
colormap parula


temp = colormap;
temp(1,:) = [1 1 1];
colormap(temp)

h = colorbar
ylabel(h,'Coherence')
caxis([0 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%