%% Demo of 2D FFAST on 'Cal' image
%
% (c) Frank Ong 2015
clc
clear
close all
setPath


%% Set Parameters
D = 2;              % number of delays (increase to get noise robustness)
noise_flag = 0;     % set noisy flag, 0 = noiseless
SNRdB = 15;         % SNR, ignore if noise_flag = 0

% Prime factors
p = [5, 7, 8];

BinsX = zeros(1,length(p));
for i = 1:length(p);
    p1 = p;
    p1(i) = 1;
    BinsX(i) = prod(p1);
end


d = length(p);

BinsY = BinsX;

nx = lcm_all(BinsX);
ny = lcm_all(BinsY);

n = nx*ny;
rho = 10^(SNRdB/10);    % adjust signal amplitude to achieve SNR
thresh = rho*0.;        % false detection threshold

%% Get Cal image
% Cal
load imCal
spectrum = zpad(imCal,[nx,ny])*sqrt(rho);
[xLoc,yLoc] = find(spectrum~=0);
k = sum(abs(spectrum(:))~=0);


%% Generate True Signal
sigma_z = sqrt(k/n)*noise_flag;
noise = sigma_z*(randn(nx,ny)+1i*randn(nx,ny))/sqrt(2);
noisy_spectrum = spectrum +  noise;
actSNR = 20*log10(norm(spectrum(:))/norm(noise(:)));
binSigma = sqrt(k./BinsX./BinsY);
binRhoSig = sqrt(rho)./binSigma;
binSNR = 20*log10(binRhoSig);

figure,imshowf(abs(noisy_spectrum),[]),title('Noisy Spectrum','FontSize',14);
signal = ifft2(noisy_spectrum);

%% Generate Bin Signals
binxSignals = zeros(D,sum(BinsX),max(BinsY));
binySignals = zeros(D,sum(BinsX),max(BinsY));
offsetX = 1;
measurements = zeros(nx,ny);
for b = 1:length(BinsX)
    binxInds = offsetX:(offsetX+BinsX(b)-1);
    binyInds = 1:BinsY(b);
    dy = 1;
    for dx = 1:D
        sampleX = dx:(nx/BinsX(b)):nx;
        sampleY = dy:(ny/BinsY(b)):ny;
        binxSignals(dx,binxInds,binyInds) = fft2(signal(sampleX,sampleY))*(nx*ny/BinsX(b)/BinsY(b));
        measurements(sampleX,sampleY) = 1;
        
    end
    dx = 1;
    for dy = 1:D
        sampleX = dx:(nx/BinsX(b)):nx;
        sampleY = dy:(ny/BinsY(b)):ny;
        binySignals(dy,binxInds,binyInds) = fft2(signal(sampleX,sampleY))*(nx*ny/BinsX(b)/BinsY(b));
        measurements(sampleX,sampleY) = 1;
        
    end
    
    offsetX = offsetX + BinsX(b);
end
m = sum(sum(measurements));


figure,imshowf(log( abs( measurements .* signal ) ), []),
title( 'Measurements','FontSize',14 );
fprintf('######## 2D FFAST #######\n');
fprintf('k = \t\t %d\n',   k);
fprintf('delays = \t %d\n',   D);
fprintf('nx = \t\t %d\n',   nx);
fprintf('ny = \t\t %d\n',   ny);
fprintf('n = \t\t %d\n',   n);
fprintf('m = \t\t %d\n',  m);
fprintf('\n');
fprintf('k/n = \t\t %.2f%%\n', k/n *100);
fprintf('m/n = \t\t %.2f%%\n',  m/n*100);
fprintf('k/m = \t\t %.2f%%\n',  k/m*100);
fprintf('\n');



for i = 1:d;
    p1 = p;
    p1(i) = 1;
    fprintf('Bin%d = \t\t ',i);
    for j = 1:d
        fprintf('%d\t',p1(j));
    end
    fprintf('x\t');
    for j = 1:d
        fprintf('%d\t',p1(j));
    end
    
    fprintf('\n');
end

if(noise_flag)
    fprintf('SNR = \t\t %.1f dB\n',SNRdB);
end
fprintf('\n');


%% Run FFAST algorithm
tic
[est_spectrum,est_xLoc,est_yLoc] = FFAST2D(binxSignals,binySignals,k,D,BinsX,BinsY,nx,ny,thresh,noise_flag);
toc

%% Show image and error

figure,imshowf(abs(est_spectrum),[0,sqrt(rho)]),title('Recovered Spectrum','FontSize',14);

Location_error = sum(sum(abs((est_spectrum~=0)-(spectrum~=0))))
NMSE = sum(sum(abs(est_spectrum-spectrum).^2))/sum(sum(abs(spectrum).^2))



