function [est_spectrum,est_locs,time] = FFAST2D(binXSignals,binYSignals,k,D,BinsX,BinsY,nX,nY,thresh,noise_flag)
% [est_spectrum,est_xLoc] = FFAST2D(binSignals,k,n,D,delays,Bins,thresh,noise_flag,kay_flag)
%
% 2D FFAST backend. Computes the sparse 2D spectrum from undersampled bin
% signals. 
%
% Inputs: 
%     binXSignals  - bin signals from undersampled DFT with x delays
%     binYSignals  - bin signals from undersampled DFT with y delays
%     k           - maximum sparsity
%     D           - number of delays per stage
%     BinsX       - x bin sizes for each stage
%     BinsY       - y bin sizes for each stage
%     nx          -  ambient x dimension
%     ny          -  ambient y dimension
%     threshold   - false detection threshold (only accept spectrum above
%     it)
%     noise_flag  - noisy or not
%
% Outputs:
%     est_spectrum     -  estimated spectrum
%     est_xLoc     -  estimated x location
%     est_xLoc     -  estimated y location
%       
% (c) Frank Ong, Sameer Pawar 2015
BinX_offset = zeros(sum(BinsX),max(BinsY));
BinX_sizes  = zeros(sum(BinsX),max(BinsY));
BinY_sizes  = zeros(sum(BinsX),max(BinsY));
T = zeros(sum(BinsX),max(BinsY));
threshold = 10^(-10);
noise_var   = (k./BinsX./BinsY)*noise_flag;
offsetX = 1;
for b = 1:length(BinsX)
    binXInds = offsetX:(offsetX+BinsX(b)-1);
    binYInds = 1:BinsY(b);
    
    T(binXInds,binYInds) = 4*(2*D-1)*noise_var(b) + threshold;
    BinX_sizes(binXInds,binYInds) = BinsX(b);
    BinY_sizes(binXInds,binYInds) = BinsY(b);
    BinX_offset(binXInds,binYInds) = offsetX-1;
    
    offsetX = offsetX + BinsX(b);
end
offsetsX = ones(1,length(BinsX));
for i = 2:1:length(BinsX)
    offsetsX(i)  =  offsetsX(i) + sum(BinsX(1:i-1));
end
time = 0;

%%
wnX = exp(1i*2*pi/nX);
wnY = exp(1i*2*pi/nY);
var_counter = 1;
est_locs = zeros(k,2);
est_vals = zeros(k,1);

tic

for it = 1:k
    singleton_found = 0;
    % Find Singletons
    signal_bins = (squeeze(sum(abs(binXSignals).^2,1))+squeeze(sum(abs(binYSignals(2:end,:,:)).^2,1))> T);

    Yx = binXSignals(:,signal_bins);
    Yy = binYSignals(:,signal_bins);
    n_signal_bins = sum(sum(signal_bins));
    [I,J] = find(signal_bins);
    signal_2_actual_mappingX = I' - BinX_offset(signal_bins)';
    signal_2_actual_mappingY = J';
    BinX_sizes_est = BinX_sizes(signal_bins)';
    BinY_sizes_est = BinY_sizes(signal_bins)';
    T_est         = T(signal_bins)';
    if (n_signal_bins)
        % ################# X ######################
        M1x = Yx(1:D-1,:,:);
        M2x = Yx(2:D,:,:);
        
        omegaX = angle(conj(M1x).*M2x);
        % find the omegas that are on boundary
        avg_positive = (median(abs(omegaX),1) > pi/2);
        omegaX  = omegaX + 2*pi*((omegaX < 0).*repmat(avg_positive,D-1,1));
        
        % apply weights
        tX  = 0:D-2;
        wtX = (1.5*D/(D^2-1))*( 1 - (  (tX-(D/2-1))/(D/2)  ).^2   );
        
        wt_omegaX = wtX*omegaX;
        wt_omegaX = wt_omegaX + 2*pi*(wt_omegaX < 0);
        est_locationX  =  mod(round(wt_omegaX*nX/(2*pi)),nX);
        % correcting the estimates using the windowing
        est_locationX = mod(round((est_locationX - (signal_2_actual_mappingX - 1))./BinX_sizes_est).*BinX_sizes_est + (signal_2_actual_mappingX - 1),nX);
        
        % ################# Y ######################
        M1y = Yy(1:D-1,:,:);
        M2y = Yy(2:D,:,:);
        
        omegaY = angle(conj(M1y).*M2y);
        % find the omegas that are on boundary
        avg_positive = (median(abs(omegaY),1) > pi/2);
        omegaY  = omegaY + 2*pi*((omegaY < 0).*repmat(avg_positive,D-1,1));
        
        % apply weights
        tY  = 0:D-2;
        wtY = (1.5*D/(D^2-1))*( 1 - (  (tY-(D/2-1))/(D/2)  ).^2   );
        
        wt_omegaY = wtY*omegaY;
        wt_omegaY = wt_omegaY + 2*pi*(wt_omegaY < 0);
        est_locationY  =  mod(round(wt_omegaY*nY/(2*pi)),nY);
        % correcting the estimates using the windowing
        est_locationY = mod(round((est_locationY - (signal_2_actual_mappingY - 1))./BinY_sizes_est).*BinY_sizes_est + (signal_2_actual_mappingY - 1),nY);
        
        
        % ########### Joint Estimate ###########
        direction_matrixX =  (wnX.^(mod((0:D-1).'*est_locationX,nX)));
        est_alphaX = sum(  Yx.*conj(direction_matrixX),1  )/D;
        direction_matrixY =  (wnY.^(mod((0:D-1).'*est_locationY,nY)));
        est_alphaY = sum(  Yy.*conj(direction_matrixY),1  )/D;
        est_alpha = (est_alphaX+est_alphaY)/2;
        est_energy = abs(est_alpha).^2;
        
        % Estimate the complete signal
        Signal_matrixX = repmat(est_alpha,[D,1,1]).*direction_matrixX;
        Signal_matrixY = repmat(est_alpha,[D,1,1]).*direction_matrixY;
        
        % Estimate noise matrix
        Noise_matrixX  =   Yx - Signal_matrixX;
        Noise_matrixY  =   Yy - Signal_matrixY;
        noise_var_est = sum(abs(Noise_matrixX).^2,1)+sum(abs(Noise_matrixY(2:end,:,:)).^2,1);
        % Set some of the singleton estimates to zero if noise is high
        singleton_est = (  noise_var_est <= T_est );   
        singleton_est = singleton_est.*(est_energy > thresh);
        
        tmp = singleton_est;
        singleton_est = zeros(sum(BinsX),max(BinsY));
        singleton_est(signal_bins) = tmp;
        tmp = est_alpha;
        est_alpha = zeros(sum(BinsX),max(BinsY));
        est_alpha(signal_bins) = tmp;
        tmp = est_locationX;
        est_locationX = zeros(sum(BinsX),max(BinsY));
        est_locationX(signal_bins) = tmp;
        tmp = est_locationY;
        est_locationY = zeros(sum(BinsX),max(BinsY));
        est_locationY(signal_bins) = tmp;
        
        % Find all singletons and remove them
        for stage = 1:1:length(BinsX)
            nx_bins = BinsX(stage);
            ny_bins = BinsY(stage);
            for binX = 1:1:nx_bins
                for binY = 1:1:ny_bins
                    bin_indexX = offsetsX(stage) + binX - 1;
                    bin_indexY = binY;
                    is_singleton = singleton_est(bin_indexX,bin_indexY);
                    if(is_singleton)
                        est_val = est_alpha(bin_indexX,bin_indexY);
                        est_locX = est_locationX(bin_indexX,bin_indexY);
                        est_locY = est_locationY(bin_indexX,bin_indexY);
                        
                        if(mod(est_locX,nx_bins) ~= (binX - 1)&&mod(est_locY,ny_bins) ~= (binY - 1))
                            continue;
                        end
                        
                        % Check if this is new singleton
                        if ( ~ismember( [est_locX,est_locY]+1, est_locs(1:var_counter-1,:),'rows' ) )
%                         if(  abs(est_spectrum(est_locX+1,est_locY+1))==0   )
                            % At least one new singleton found
                            singleton_found = 1;
                            if(var_counter <= k)
                                est_locs(var_counter,1) = est_locX + 1;
                                est_locs(var_counter,2) = est_locY + 1;
                                est_vals(var_counter) = est_val;
%                                 est_spectrum(est_locX+1,est_locY+1) = est_val;
                                
                                var_counter = var_counter + 1;
                            end
                            
                            loc_signalX = est_val.*(wnX.^(mod((0:D-1).'*est_locX,nX)));
                            loc_signalY = est_val.*(wnY.^(mod((0:D-1).'*est_locY,nY)));
                            for b = 1:length(BinsX)
                                % Find the bin locations where it hashes
                                hashX = mod(est_locX,BinsX(b))+offsetsX(b);
                                hashY = mod(est_locY,BinsY(b))+1;
                                binXSignals(:,hashX,hashY) = binXSignals(:,hashX,hashY) - loc_signalX;
                                binYSignals(:,hashX,hashY) = binYSignals(:,hashX,hashY) - loc_signalY;
                            end
                        end
                    end
                end
            end
            
            if(var_counter >= k)
                break;
            end
        end
        
    else
        break;
    end
    
    if(~singleton_found)
        break;
    end
    
end

time = toc;

% [est_xLoc,est_yLoc] = find(est_spectrum~=0);
% est_xLoc = est_xLoc';
% est_yLoc = est_yLoc';


est_spectrum = zeros(nX,nY);
for i = 1:var_counter-1
    est_spectrum(est_locs(i,1), est_locs(i,2)) = est_vals(i);
end
