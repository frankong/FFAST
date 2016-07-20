%% Demo of 2D FFAST on the 'House' image using the good thomas mapping
%
% (c) Frank Ong 2015
clc
clear
close all
setPath

%% Set Parameters
D = 3;      % number of delays for each chain
C = 1;      % number of chains for delays
SNRdB = 25;
noise_flag = 1;
kay_flag = 1; % Use kay's estimation
very_sparse = 0;

N = D * C;  % total number of delays

if (~kay_flag)
    D = N;
    C = 1;
end

rho = 10^(SNRdB/10);


thresh = rho * 0.5;

p = [13, 19, 14, 17];

if (very_sparse)
    Bins = p;
else
    Bins = zeros(1,length(p));
    for i = 1:length(p);
        p1 = p;
        p1(i) = 1;
        Bins(i) = prod(p1);
    end
end

n = lcm_all(Bins);
nx = prod(p(1:ceil(end/2)));
ny = prod(p(ceil(end/2)+1:end));


%% Generate Delays

rand_loc = 1;

if (kay_flag)
    delays = zeros(1,N);
    r = 1;
    for i = 1:C
        if rand_loc
            r = randi(n)-1;
        end
        skip = 2^(i-1);
        delays((i-1)*D+(1:D)) = r + (0:skip:(skip*D-1));
        if ~rand_loc
            r = r + i*D;
        end
    end
    
else
    delays = randi(n,1,N);
end


%% Get 2D Spectrum

load imHouse
spectrum2 = crop(im,nx,ny);
[xLoc,yLoc] = find(spectrum2~=0);
k = sum(abs(spectrum2(:))~=0);

scale = 1 / max(abs(spectrum2(:))) * sqrt(rho);

spectrum2 = spectrum2 * scale;


figure,imshowf(abs(spectrum2),[0,sqrt(rho)]),title('Original Spectrum','FontSize',14),drawnow;


%% Generate 2D Signal

sigma_z = sqrt(k/n)*noise_flag;
noise = sigma_z*(randn([nx,ny])+1i*randn([nx,ny]))/sqrt(2);

noisy_spectrum2 = spectrum2 + noise;

signal2 = ifft2(noisy_spectrum2);


figure,imshowf(abs(noisy_spectrum2),[0,sqrt(rho)]),title('Noisy Spectrum','FontSize',14),drawnow;


%% Generate 1D signal
signal = crtmapping(signal2);


%% Generate Bin Signals
binSignals = zeros(N,sum(Bins));
offset = 1;
measurements = zeros(1,n);
for f = 1:length(Bins)
    binInds = offset:(offset+Bins(f)-1);
    for dind = 1:N
        shift_signal = circshift( signal, [0,-delays(dind)] );
        sub_signal = shift_signal(1:(n/Bins(f)):end);
        
        binSignals(dind,binInds) = fft(sub_signal)*(n/Bins(f));
        
        measurements = circshift( measurements, [0,-delays(dind)] );
        measurements(1:(n/Bins(f)):end) = 1;
        measurements = circshift( measurements, [0,delays(dind)] );
    end
    offset = offset + Bins(f);
end
m = sum(measurements);
measurements2 = goodmapping(measurements,[nx,ny]);
calcSNR = 10*log10(2*n^2*log(n)/(N^2*(m/1.23)^2));

fprintf('######## 2D FFAST #######\n');
fprintf('k = \t\t %d\n',   k);
fprintf('D = \t\t %d\n',   D);
fprintf('C = \t\t %d\n',   C);
fprintf('N = \t\t %d\n',   N);
fprintf('n = \t\t %d\n',   n);
fprintf('m = \t\t %d\n',  m);
fprintf('\n');
fprintf('k/n = \t\t %.2f%%\n', k/n *100);
fprintf('m/n = \t\t %.2f%%\n',  m/n*100);
fprintf('k/m = \t\t %.2f%%\n',  k/m*100);
fprintf('\n');


if (very_sparse)
    for i = 1:length(p);
        fprintf('Bin%d = \t\t ',i);
        fprintf('%d\t',p(i));
        fprintf('\n');
    end
    fprintf('\n');
    
else
    for i = 1:length(p);
        p1 = p;
        p1(i) = 1;
        fprintf('Bin%d = \t\t ',i);
        for j = 1:length(p)
            fprintf('%d\t',p1(j));
        end
        
        fprintf('\n');
    end
    fprintf('\n');
end

for i = 1:C
    
    fprintf('Delay%d = \t ',i);
    for j = 1:D
        fprintf('%d\t',delays( j + (i-1)*D ));
    end
    
    fprintf('\n');
end
fprintf('\n');


if(noise_flag)
    fprintf('SNR = \t\t %.1f dB\n',SNRdB);
end
fprintf('\n');

figure,imshowf( log( abs( measurements2 .* signal2 ) ), [] ),
title('Measured Signal','FontSize',14);


%% Estimate Spectrum
tic
[est_spectrum,est_xLoc] = FFAST1D(binSignals,k,n,D,delays,Bins,thresh,noise_flag,kay_flag);
toc

est_spectrum2 = goodmapping(est_spectrum,[nx,ny]);

Location_error = sum(sum(abs((est_spectrum2~=0)-(spectrum2~=0))))
NMSE = sum(sum(abs(est_spectrum2-spectrum2).^2))/sum(sum(abs(spectrum2).^2))

figure,imshowf(abs(est_spectrum2),[0,sqrt(rho)*0.1]),title('Recovered Spectrum','FontSize',14);

