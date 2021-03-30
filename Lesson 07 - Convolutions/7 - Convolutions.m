% Lesson 07 - Convolutions 

% Short demonstration of a convolution

X = [0 2 1 1 0]
K = [1 2 2]

ConvXK = conv(X,K)

%% Restricting convolution product size to the same as the signal

X = [0 2 1 1 -1 7 -2 -3 8 1 -1 2 -3 4]
K = [1 1 -2]

ConvXK = conv(X,K)
ConvXKB = conv(X,K,'same')

subplot(311)
plot(2:15,X,'bo-')
xlim([0 17])

subplot(312)
plot(K,'ro-')
xlim([0 17])

subplot(313)
plot(ConvXK,'ko-')
hold on
plot(2:15,ConvXKB,'yo-')
hold off
xlim([0 17])

%% Example application of convolutions

clf

srate = 1000;
dt = 1/srate;
t = dt:dt:2;

LFP = randn(size(t));
subplot(211)
plot(t,LFP)

order = 500
% moving average 
K = ones(1,order)/order;

% Convola = conv(LFP,K);
Convol = conv(LFP,K,'same');

hold on
plot(t,Convol,'k-','linew',3)
hold off

subplot(212)
[Pxx F] = pwelch(Convol,length(Convol),[],2^16,srate);
hold off
plot(F,Pxx)
xlim([0 500])
xlabel('freq (hz)')



%% Convolution effects on the frequency domain

subplot(312)

K = ones(1,50)

[Pkk F] = pwelch(K,length(K),[],2^16,srate);
hold off
plot(F,Pkk)
xlim([0 500])
xlabel('freq (hz)')
ylabel('Pkk')
% %

subplot(311)
[Pxx F] = pwelch(LFP,length(K),[],2^16,srate);
hold off
plot(F,Pxx)
xlim([0 500])
xlabel('freq (hz)')
ylabel('Pxx')


%% Varying kernel size

% clf

srate = 1000;
dt = 1/srate;
t = dt:dt:2;

% LFP = randn(size(t));
% LFP(500:1500) = LFP(500:1500) + 30*exp(-t(500:1500));

order = 50;

K = ones(1,order)/order;

Convol = conv(LFP,K,'same')

subplot(111)
plot(t,LFP)
hold on
plot(t,Convol,'r-','linew',3)
plot([0.2 0.2],[-5 20],'r--','linew',2)
hold off

xlabel('time (s)')

%% Again

clf
clear
clc

% %
srate = 1000;
dt = 1/srate;
t = dt:dt:1;

LFP= sin(2*pi*10*t)+2*sin(2*pi*0.5*t)+0.4*sin(2*pi*40*t);
LFP = LFP + randn(size(t));

order = 200;
K = sin(2*pi*10*t(1:order));
norm = sum(K.^2); % to avoid adding or removing energy

% %

figure(2)

subplot(311)
[Pxx F] = pwelch(LFP,length(LFP),[],2^16,srate);
hold off
plot(F,Pxx)
xlim([0 50])
xlabel('freq (hz)')
ylabel('Pxx')

subplot(312)
[Pkk F] = pwelch(K/norm,length(K),[],2^16,srate);
hold off
plot(F,Pkk)
xlim([0 50])
xlabel('freq (hz)')
ylabel('Pkk')

subplot(313)
plot(F,Pkk.*Pxx)
xlim([0 50])
xlabel('freq (hz)')
ylabel('Power convol')


%% Zero-padding

figure(1)
% padding with zeros
LFP = [zeros(1,order) LFP zeros(1,order)];
t = (1:length(LFP))*dt;

for j = 0:1:length(LFP)-order
%     j = 200    
subplot(211)
plot(t,LFP)
hold on
plot(t((1:order)+j),K(order:-1:1),'r-','linew',3)
hold off
ylim([-4 4])

subplot(212)

plot(t(order/2+j),sum(K(order:-1:1).*LFP((1:order)+j))/norm,'ko')
hold on
xlim([0 t(end)])
ylim([- 4 4])
pause(0.001)
end

%% Dealing with phase shifts

figure(1)
Convol = conv(LFP,K/norm,'same')

% to avoid phase shifts, convolve the 
% convolution with the inverted kernel
Convol2 = conv(Convol,K(end:-1:1)/norm,'same')

subplot(212)
hold on
plot(t,Convol,'r-')
plot(t,Convol2,'g-')
hold off

subplot(211)
hold on
plot(t,Convol,'r-')
plot(t,Convol2,'g-')
hold off



%%
figure(2)
subplot(313)

[Pkk F] = pwelch(Convol2,length(Convol),[],2^16,srate);
hold on
plot(F,Pkk)
xlim([0 50])
xlabel('freq (hz)')
ylabel('Pkk')






































