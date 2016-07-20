%% Demo of 2D FFAST on synthetic sparse images
%
% (c) Frank Ong 2015
clc
clear
close all
setPath

niter = 1000;

k_array = 100*(1:3); % sparsity
scale_array = 1:2:20; % scaling for n

size_array = zeros(length(scale_array), length(k_array));
time_array = zeros(length(scale_array), length(k_array));

D = 2;              % number of delays (increase to get noise robustness)
noise_flag = 0;     % set noisy flag
SNRdB = 15;         % SNR, ignore if noise_flag = 0
rho = 10^(SNRdB/10);    % adjust signal amplitude to achieve SNR
thresh = rho*0.99;        % false detection threshold

%%

for k_ind = 1:length(k_array)
    
    k = k_array(k_ind);
    
    for s_ind = 1:length(scale_array)
        
        
        %% Set Parameters
        
        % Prime factors
        p = [5,7,9];
        d = length(p);
        
        BinsX = zeros(1,length(p));
        for i = 1:length(p);
            p1 = p;
            p1(i) = 1;
            BinsX(i) = prod(p1);
        end
        
        BinsY = BinsX;
        
        
        
        nx = lcm_all(BinsX) * scale_array(s_ind);
        ny = lcm_all(BinsY);
        
        n = nx*ny;
        size_array(s_ind, k_ind) = n;
        
        %% Generate True Spectrum
        % Random
        spectrum = zeros(nx,ny);
        loc = randperm(n,k);
        xLoc = mod(loc,nx)+1;
        yLoc = ceil(loc/nx);
        spectrum(xLoc+(yLoc-1)*nx) = sqrt(rho)*exp(1i*rand(1,k)*2*pi);
        
        
        %% Generate True Signal
        sigma_z = sqrt(k/n)*noise_flag;
        noise = sigma_z*(randn(nx,ny)+1i*randn(nx,ny))/sqrt(2);
        noisy_spectrum = spectrum +  noise;
        actSNR = 20*log10(norm(spectrum(:))/norm(noise(:)));
        binSigma = sqrt(k./BinsX./BinsY);
        binRhoSig = sqrt(rho)./binSigma;
        binSNR = 20*log10(binRhoSig);
        
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
        
        
        
        fprintf('######## 2D FFAST #######\n');
        fprintf('k = \t\t %d\n',   k);
        fprintf('nx = \t\t %d\n',   nx);
        fprintf('ny = \t\t %d\n',   ny);
        fprintf('n = \t\t %d\n',   n);
        fprintf('m = \t\t %d\n',  m);
        fprintf('\n');
        
        
        
        %% Estimate Spectrum
        for it = 1:niter
            [est_spectrum,est_locs,time] = FFAST2D(binxSignals,binySignals,k,D,BinsX,BinsY,nx,ny,thresh,noise_flag);
            time_array(s_ind, k_ind) = time_array(s_ind, k_ind) + time;
%             NMSE = sum(sum(abs(est_spectrum-spectrum).^2))/sum(sum(abs(spectrum).^2))
        end
        
        
    end
    
end

time_array = time_array / niter;

%%
figure,plot(size_array/315,time_array,'LineWidth',5)
legend('k = 100', 'k = 200', 'k = 300');
xlabel('Signal Dimension Nx','FontSize',14)
ylabel('Average Run Time (sec)','FontSize',14)
set(gca,'FontSize',14)

