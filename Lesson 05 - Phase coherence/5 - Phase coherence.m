% Lesson 5 - Coherence

% Understanding phase coherence

clear
R = 3;
theta = deg2rad(45);
Fx = R*exp(1i*theta);

R = 0.5;
theta = deg2rad(145);
Fy = R*exp(1i*theta);

normFx = Fx/abs(Fx);
normFy = Fy/abs(Fy);

plot([real(Fx)],[imag(Fx)],'bo','markerf','b')
hold on
plot([0 real(Fx)],[0 imag(Fx)],'b-')
plot([real(normFx)],[imag(normFx)],'bo','markerf','b')
plot([0 real(normFx)],[0 imag(normFx)],'b-')

plot([real(Fy)],[imag(Fy)],'ro','markerf','r')
hold on
plot([0 real(Fy)],[0 imag(Fy)],'r-')
plot([real(normFy)],[imag(normFy)],'ro','markerf','r')
plot([0 real(normFy)],[0 imag(normFy)],'r-')

plot(exp(1i*(0:0.01:2*pi)),'k-')
hold off

text([real(Fx)+0.1],[imag(Fx)],'Fx')
text([real(Fy)+0.1],[imag(Fy)],'Fy')

text([real(normFx)+0.1],[imag(normFx)],'Fx/|Fx|')
text([real(normFy)+0.1],[imag(normFy)],'Fy/|Fy|')

xlim([-3 3])
ylim([-3 3])

axis square

%% Averaging unitary vectors

clear Vall
for nvector = 1:20
    theta = pi/2 + 1*randn;
    v = exp(i*theta);
    
    Vall(nvector) = v;
    
    plot([real(v)],[imag(v)],'bo','markerf','b')
    hold on
    plot([0 real(v)],[0 imag(v)],'b-')
    plot(exp(1i*(0:0.01:2*pi)),'k-')
    
    % hold off
    
    xlim([-1 1])
    ylim([-1 1])
    
    axis square
end

MeanVector = mean(Vall);
plot([real(MeanVector)],[imag(MeanVector)],'ko','markerf','k')
hold on
plot([0 real(MeanVector)],[0 imag(MeanVector)],'k-')
hold off

% coherence is defined as the length of
% the mean vector (mean over all unitary vectors)
coherence = abs(MeanVector);

title(['Coherence = ' num2str(coherence)])


%% From the beggining ...

clear
clf

srate = 1000;
dt = 1/srate;
t = dt:dt:1;

phi = -deg2rad(45);

X = 3*sin(2*pi*8*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*8*t+phi)+0.3*randn(size(t));

subplot(311)
plot(t,X)
hold on
plot(t,Y)
hold off
xlim([0 1])

subplot(3,1,[2 3])

f = 8; % Frequency of Fourier Kernel
K = exp(-1i*2*pi*f*t);

FX = mean(X.*K); % actually, FX is the sum, we are using the mean just for visualization
FY = mean(Y.*K);

nFX = FX/abs(FX);
nFY = FY/abs(FY);

nFXY = nFX*conj(nFY);
% nFXY = FX*conj(FY)/(abs(FX)*abs(FY));

plot([real(FX)],[imag(FX)],'bo')
hold on
plot([0 real(FX)],[0 imag(FX)],'b-')
plot([real(FY)],[imag(FY)],'ro')
plot([0 real(FY)],[0 imag(FY)],'r-')

plot([real(nFX)],[imag(nFX)],'bo','markerf','b')
plot([0 real(nFX)],[0 imag(nFX)],'b-')
plot([real(nFY)],[imag(nFY)],'ro','markerf','r')
plot([0 real(nFY)],[0 imag(nFY)],'r-')

plot([real(nFXY)],[imag(nFXY)],'ko','markerf','k')
plot([0 real(nFXY)],[0 imag(nFXY)],'k-')

plot(exp(1i*(0:0.01:2*pi)),'k-')
hold off

xlim([-2 2])
ylim([-2 2])
axis square

%% Loop across time windows

clear
clf

srate = 1000;
dt = 1/srate;
t = dt:dt:100;

f = 7.5; % Frequency of Fourier Kernel
K = exp(-1i*2*pi*f*t);

clear nFXYAll

for nwindow = 1:1
    
    phi = -deg2rad(0)+0.3*randn;
    X = 3*sin(2*pi*8*t)+0.3*randn(size(t));
    Y = 0.5*sin(2*pi*8*t+phi)+0.3*randn(size(t));
    
    FX = mean(X.*K); % actually, FX is the sum, we are using the mean just for visualization
    FY = mean(Y.*K);
    
    nFX = FX/abs(FX);
    nFY = FY/abs(FY);
    
    nFXY = nFX*conj(nFY);
    
    nFXYAll(nwindow) = nFXY;
    % nFXY = FX*conj(FY)/(abs(FX)*abs(FY));
    
    % % plot([real(FX)],[imag(FX)],'bo')
    % % hold on
    % % plot([0 real(FX)],[0 imag(FX)],'b-')
    % % plot([real(FY)],[imag(FY)],'ro')
    % % plot([0 real(FY)],[0 imag(FY)],'r-')
    % %
    % % plot([real(nFX)],[imag(nFX)],'bo','markerf','b')
    % % plot([0 real(nFX)],[0 imag(nFX)],'b-')
    % % plot([real(nFY)],[imag(nFY)],'ro','markerf','r')
    % % plot([0 real(nFY)],[0 imag(nFY)],'r-')
    
    
    plot([real(nFXY)],[imag(nFXY)],'ko','markerf','k')
    hold on
    plot([0 real(nFXY)],[0 imag(nFXY)],'k-')
    
    plot(exp(1i*(0:0.01:2*pi)),'k-')
    hold on
    
    xlim([-2 2])
    ylim([-2 2])
    axis square
    
    title(['Nwindow = ' num2str(nwindow)])
    pause(0.005)
end


Cxy = abs(mean(nFXYAll));

text(0,1.5,['Coherence = ' num2str(Cxy)],'fontsize',16)

hold on

plot([real(mean(nFXYAll))],[imag(mean(nFXYAll))],...
    'ko','markerf','y')
hold on
plot([0 real(mean(nFXYAll))],...
    [0 imag(mean(nFXYAll))],'y-','linew',3)

hold off

%% Coherence spectrum

srate = 1000;
dt = 1/srate;
t = dt:dt:1;

clear CxySpectrum

freqvector = 0:0.1:10;
count = 0;

for f = freqvector
    count = count+1;
    
    K = exp(-1i*2*pi*f*t);
    
    clear nFXYAll
    
    for nwindow = 1:100
        
        phi = -deg2rad(90)+0*randn;
        X = 3*sin(2*pi*5*t)+0.3*randn(size(t));
        Y = 0.5*sin(2*pi*5*t+phi)+0.3*randn(size(t));
        
        % FX = mean(X.*K); % actually, FX is the sum, we are using the mean just for visualization
        % FY = mean(Y.*K);
        
        W = hamming(length(X));
        FX = mean((W.*X).*K); % actually, FX is the sum, we are using the mean just for visualization
        FY = mean((W.*Y).*K);
        
        nFX = FX/abs(FX);
        nFY = FY/abs(FY);
        
        nFXY = nFX*conj(nFY);
        
        nFXYAll(nwindow) = nFXY;
        
        plot([real(nFXY)],[imag(nFXY)],'ko','markerf','k')
        hold on
        plot([0 real(nFXY)],[0 imag(nFXY)],'k-')
        plot(exp(1i*(0:0.01:2*pi)),'k-')
        hold on
        
        xlim([-1 1])
        ylim([-1 1])
        axis square
        
        title(['Nwindow = ' num2str(nwindow)])
        % pause(0.005)
    end
    
    Cxy = abs(mean(nFXYAll));
    
    CxySpectrum(count) = Cxy;
    
    title(['Coherence = ' num2str(Cxy) '  Frequency = ' num2str(f) ' Hz'])
    
    plot([real(mean(nFXYAll))],[imag(mean(nFXYAll))],...
        'ko','markerf','y')
    hold on
    plot([0 real(mean(nFXYAll))],...
        [0 imag(mean(nFXYAll))],'y-','linew',3)
    
    % pause
    
    hold off
end


%%

plot(freqvector,CxySpectrum)
hold on
plot(freqvector,CxySpectrum.^2)
hold off

xlabel('Frequency (Hz)')
ylabel('Coherence')


%% Using Matlab built-in function

% mscohere % ms = magnitude squared

srate = 1000;
dt = 1/srate;
t = dt:dt:10;

phi = -deg2rad(90)+0.5*randn;
X = 3*sin(2*pi*10*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*10*t+phi)+0.3*randn(size(t));

windowlength = 2*srate;
overlap = 0;
freqvector = 0:0.1:20;

% [Cxy F] = mscohere(X,Y,windowlength,overlap,freqvector,srate)
% plot(F,Cxy)
% xlabel('Frequency (Hz)')
% ylabel('Coherence')

nfft = 2^16;
[Cxy F] = mscohere(X,Y,windowlength,overlap,nfft,srate);
hold on
plot(F,Cxy)
xlabel('Frequency (Hz)')
ylabel('Coherence')

xlim([0 20])

%% Coherogram

srate = 1000;
dt = 1/srate;
t = dt:dt:200;

phi = -deg2rad(90)+0.5*randn;
X = 3*sin(2*pi*10*t)+0.3*randn(size(t));
Y = 0.5*sin(2*pi*10*t+phi)+0.3*randn(size(t));

win = 10*srate; % window for Cxy computation
step = 0.1*win;
Nwin = (length(t)-win)/step+1;

cohwin = 2*srate;
nfft = 2^16;
overlap = 0;

for nwin = 1:Nwin
    winidx = (1:win) + (nwin-1)*step;
    [Cxy F] = mscohere(X(winidx),Y(winidx),cohwin,overlap,nfft,srate);
    Coherogram(nwin,:) = Cxy;
    T(nwin) = t(winidx(win/2));
end

%%

clf
imagesc(T,F,Coherogram')
axis xy
ylim([0 20])
xlabel('Time (s)')
ylabel('Frequency (Hz)')

h = colorbar;
h.Label.String = 'Coherence';









