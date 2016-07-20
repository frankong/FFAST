%% Demo of 2D FFAST on 'Brain' image
%
% (c) Frank Ong 2015
clc
clear
close all
setPath


%% Set Parameters
D = 15;
noise_flag = 1;
kay_flag = 0;
SNRdB = 20;

p = [7,8,9];

l = length(p);
BinsX = zeros(1,length(p));
for i = 1:length(p);
    p1 = p;
    p1(i) = 1;
    BinsX(i) = prod(p1);
end
BinsY = BinsX;

nx = lcm_all(BinsX);
ny = lcm_all(BinsY);

n = nx*ny;
rho = 10^(SNRdB/10);
thresh = rho*0.02;

%% Get Actual Spectrum
load brain7t
im = crop(im,nx,ny);
linPhase = repmat(1-exp(1i*2*pi/nx*(0:nx-1)'),[1,ny])/sqrt(2);
spectrum = fft2(ifft2(im).*linPhase);

scale = sqrt(rho)/max(abs(spectrum(:)));
spectrum = spectrum*scale;
[xLoc,yLoc] = find(abs(spectrum)<=sqrt(thresh));
k = sum(sum(abs(spectrum)>=sqrt(thresh)));

figure,imshowf(abs(im),[]),title('Original image','FontSize',14);
figure,imshowf(abs(spectrum),[]),title('Original edge image','FontSize',14);

%% Get Actual Signal
signal = ifft2(spectrum);

%% Get Random Delay
delay = zeros(D,2,length(BinsX));

for stage = 1:length(BinsX)
    
    r = randperm( nx*ny/BinsX(stage)/BinsY(stage), D );
    x = mod( r-1 , nx/BinsX(stage) )';
    y = floor( (r-1) / (nx/BinsX(stage)) )';
    
    delay(:,1,stage) = x;
    delay(:,2,stage) = y;
    delay(:,:,stage) = sortrows(delay(:,:,stage));
    
end

%% Generate Bin Signals
binSignals = zeros(D,sum(BinsX),max(BinsY));
offsetX = 1;
measurements = zeros(nx,ny);

% Sample
for stage = 1:length(BinsX)
    binxInds = offsetX:(offsetX+BinsX(stage)-1);
    binyInds = 1:BinsY(stage);
    
    for d = 1:D
        dx = delay(d,1,stage);
        dy = delay(d,2,stage);
        sampleX = (dx+1):(nx/BinsX(stage)):nx;
        sampleY = (dy+1):(ny/BinsY(stage)):ny;
        binSignals(d,binxInds,binyInds) = fft2(signal(sampleX,sampleY))*(nx*ny/BinsX(stage)/BinsY(stage));
        measurements(sampleX,sampleY) = 1;
        
    end
    
    offsetX = offsetX + BinsX(stage);
end

% Get fully sampled low frequency signal
lowSize = 20;
sig = ifft2( im ) * scale;
lowFreqRegion = zeros(nx,ny);
lowFreqRegion(1,:) = 1;
lowFreqRegion(1:lowSize,:) = 1;
lowFreqRegion(end-lowSize+1:end,:) = 1;
lowFreqRegion = lowFreqRegion>0;
lowFreqSignal = sig( lowFreqRegion );

measurements(lowFreqRegion) = 1;

figure,imshowf(fftshift(measurements));
title('Measurements','FontSize',14);

m = sum(sum(measurements));


fprintf('######## 2D RFFAST #######\n');
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



for i = 1:l
    p1 = p;
    p1(i) = 1;
    fprintf('Bin%d = \t\t ',i);
    for j = 1:l
        fprintf('%d\t',p1(j));
    end
    fprintf('x\t');
    for j = 1:l
        fprintf('%d\t',p1(j));
    end
    
    fprintf('\n');
end

fprintf('\n');


%% Estimate Spectrum
[est_spectrum, num_it] = RFFAST2D(k,noise_flag,kay_flag,binSignals,delay,BinsX,BinsY,nx,ny,thresh);

%% Show results and errors

figure,imshowf(abs(est_spectrum),[]),title('Recovered Edge Spectrum','FontSize',14);

NMSE = sum(sum(abs(est_spectrum-spectrum).^2))/sum(sum(abs(spectrum).^2))


%% Invert edge transform to obtain the full 'Brain' image

full_sig = ifft2( est_spectrum );
full_sig = full_sig./linPhase;
full_sig(lowFreqRegion) = lowFreqSignal;

est_full_spectrum = fft2( full_sig );


figure,imshowf(abs(est_full_spectrum),[]);
,title('Recovered Full Spectrum','FontSize',14);

