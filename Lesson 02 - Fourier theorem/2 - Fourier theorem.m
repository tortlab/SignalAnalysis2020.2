% Complex number representation in the complex plane

clf
clc

x = 3;
y = 2;

z = x + 1i*y

% Polar representation

R = sqrt(x^2 + y^2);
theta = atan(y/x);

z2 = R*exp(1i*theta)

plot([0 real(z2)],[0 imag(z2)],'k-')
hold on
plot(real(z2),imag(z2),'ko','markerf','y','markers',14)
hold off

xlim([-5 5])
ylim([-5 5])

axis square

%% Again - different notation

R = 4
theta = deg2rad(-45)
z = R*exp(1i*theta)


plot([0, real(z)],[0, imag(z)],'k-')
hold on
plot(real(z),imag(z),'ko','markerf','y','markers',14)
hold off

xlim([-5 5])
ylim([-5 5])

axis square


%% Animating the vector

for theta = 0:0.01:2*pi

R = 4;
z = R*exp(-1i*theta);

plot([0, real(z)],[0, imag(z)],'k-')
hold on
plot(real(z),imag(z),'ko','markerf','y','markers',14)
hold off

xlim([-5 5])
ylim([-5 5])

axis square

pause(0.001)

end

%% 

% Fourier kernel:
% exp(-i*2*pi*f*t)

f = 10

for t=0:0.01:1

theta = 2*pi*f*t;
R = 4;
z = R*exp(-1i*theta);

plot([0, real(z)],[0, imag(z)],'k-')
hold on
plot(real(z),imag(z),'ko','markerf','y','markers',14)
hold off

xlim([-5 5])
ylim([-5 5])

axis square

title(['Time = ' num2str(t) ' s'])

pause(0.01)

end

%% Fake LFP (just a senoid)

f = 5
tvector = 0:0.001:1;
LFP = sin(2*pi*f*tvector);


LFP = sin(2*pi*3*tvector)+0.2*sin(2*pi*14*tvector);

plot(tvector,LFP)

xlabel('Time (s)')
ylabel('mV')
set(gcf,'color','w')

%% Fourier transform

clf
% sum LFP(t)*exp(-i*2*pi*f*t)

Tmax = 1;

% f1 = 3 % signal frequency,in Hz
f2 = 0 % frequency of the Fourier Kernell, in Hz

Zall = [] % this will accumulate all products

for t=0:0.01:Tmax

% LFP = abs(sin(2*pi*f1*t));

% LFP = (sin(2*pi*f1*t+pi));

% LFP = 2*sin(2*pi*2*t)+sin(2*pi*10*t+pi);

LFP = 2+sin(2*pi*2*t);

subplot(3,1,1)
plot(t,LFP,'k.')
hold on
xlim([0 Tmax])
ylim([-3 3])

    title(['Time = ' num2str(t) ' s'])
% % % 
subplot(3,1,[2 3])
    
    % Fourier Kernel
    %%%%%%%%%%%%%%%%%%%
    z = exp(-1i*2*pi*f2*t);
    %%%%%%%%%%%%%%%%%%%%
    
    % Signal x fourier kernel
    
    Z = LFP*z;
    
    Zall = [Zall,Z];
    
    plot([0, real(Z)],[0, imag(Z)],'k-')
    hold on
    plot(real(Z),imag(Z),'ko','markerf','y','markers',7)
% hold off


    % plotting the unitary circle
    plot(exp(1i*(0:0.001:2*pi)))
    
    hold on
%     plot([0 real(z)],[0 imag(z)],'b-')
%     hold off
    
    xlim([-3 3])
    ylim([-3 3])
    axis square
    
    pause(0.01)
end
 
% %
% computing the fourier transform
% FX = sum(Zall); % this is the correct definition

FX = mean(Zall); % using the mean just for visualization


% Power = abs(FX)^2
Power = FX*conj(FX)

plot([0, real(FX)],[0, imag(FX)],'r-')
    hold on
    plot(real(FX),imag(FX),'ro',...
        'markerf','r','markers',7)

    title(['Power at frequency ' num2str(f2) ...
       ' Hz = ' num2str(Power) ] )
    
subplot(3,1,1)
hold off

%% Power spectrum

clf

Tmax = 1;

% f1 = 3 % signal frequency,in Hz
% f2 = 0 % frequency of the Fourier Kernell, in Hz

clear PowerSpectrum

FreqVector = 1:1:10;

for f2 = FreqVector

Zall = []; % this will accumulate all products

for t=0:0.01:Tmax

LFP = sin(2*pi*9*t);

if t>0.5
    LFP = 0;
end


subplot(3,1,1)
plot(t,LFP,'k.')
hold on
xlim([0 Tmax])
ylim([-3 3])

    title(['Time = ' num2str(t) ' s'])
% % % 
subplot(3,1,[2 3])
    
    % Fourier Kernell
    %%%%%%%%%%%%%%%%%%%
    z = exp(-1i*2*pi*f2*t);
    %%%%%%%%%%%%%%%%%%%%
    
    % Signal x fourier kernel
    
    Z = LFP*z;
    
    Zall = [Zall,Z];
    
    plot([0, real(Z)],[0, imag(Z)],'k-')
    hold on
    plot(real(Z),imag(Z),'ko','markerf','y','markers',7)
% hold off


    % plotting the unitary circle
    plot(exp(1i*(0:0.001:2*pi)))
    
    hold on
%     plot([0 real(z)],[0 imag(z)],'b-')
%     hold off
    
    xlim([-3 3])
    ylim([-3 3])
    axis square
    
    pause(0.001)
end
 
% %
% computing the fourier transform
% FX = sum(Zall); % this is the correct definition

FX = mean(Zall); % using the mean just for visualization


% Power = abs(FX)^2
Power = FX*conj(FX);

plot([0, real(FX)],[0, imag(FX)],'r-')
    hold on
    plot(real(FX),imag(FX),'ro',...
        'markerf','r','markers',7)

    title(['Power at frequency ' num2str(f2) ...
       ' Hz = ' num2str(Power) ] )
    
   pause(0.5)
 PowerSpectrum(f2)=Power;
 hold off
 
end
   
%%

subplot(3,1,[2 3])
plot(FreqVector,PowerSpectrum,'ko-','linewidth',3,...
    'markersize',14,'markerfacecolor','w')


stem(FreqVector,PowerSpectrum,'k-','linew',3)

xlabel('Frequency (Hz)')
ylabel('Power')
set(gca,'fontsize',14)



% % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% a little detour about plotting


% plotting a dot

x = 3
y = 2

plot(x,y,'ro','markerfacecolor','y',...
    'markersize',20)

% matlab uses RGB system for colors

plot(x,y,'ro','markerfacecolor',[0 0 0],...
    'markersize',20)

%%
% plotting more than one point

% plot([1 2 3 5],[4 2 -1 2])

plot([1 2 3 5],[4 2 -1 2],'ysq-', ...
    'markersize',20,'markerfacecolor','m')
hold on
plot(1:5,[1:5].^2)
hold off

axis square

% xlim([0 6])
% ylim([-2 6])


%% Working with subplots

clf

% subplot(Rows,Columns,PanelNumber)

subplot(2,2,1)
plot([1 2 3 5],[4 2 -1 2],'ysq-', ...
    'markersize',20,'markerfacecolor','m')


subplot(2,2,2)
plot(1:5,[1:5].^2)


subplot(2,1,2)
plot(1:5,sqrt([1:5]))

% 
% subplot(2,2,4)
% plot(1:5,exp([1:5]))