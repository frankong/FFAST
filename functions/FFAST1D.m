function [est_spectrum,est_xLoc] = sFFT1D(binSignals,k,n,D,delays,Bins,thresh,noise_flag,kay_flag)
% [est_spectrum,est_xLoc] = sFFT1D(binSignals,k,n,D,delays,Bins,thresh,noise_flag,kay_flag)
%
% 1D FFAST backend. Computes the sparse 1D spectrum from undersampled bin
% signals. 
%
% Inputs: 
%     binSignals  - bin signals from undersampled DFT
%     k           - maximum sparsity
%     n           - ambient dimension
%     D           - number of delays per stage
%     delays      - delay numbers for each stage
%     Bins        - bin sizes for each stage
%     threshold   - false detection threshold (only accept spectrum above
%     it)
%     noise_flag  - noisy or not
%     kay_flag    - use kay estimator or not 
%
% Outputs:
%     est_spectrum     -  estimated spectrum
%     est_xLoc     -  estimated location
%       
% (c) Frank Ong, Sameer Pawar 2015


N = length(delays);
B = N / D;
% n = lcm_all(Bins);

%% Get Kay weights
t  = 0:D-2;
wt = (1.5*D/(D^2-1))*( 1 - (  (t-(D/2-1))/(D/2)  ).^2   );


%% Initialize indicies
Bin_offset = zeros(1,sum(Bins));
Bin_sizes  = zeros(1,sum(Bins));
T = zeros(1,sum(Bins));
threshold = 10^(-10);
noise_var   = (k./Bins)*noise_flag;
offset = 1;
for b = 1:length(Bins)
    binInds = offset:(offset+Bins(b)-1);
    T(binInds) = 4*N*noise_var(b) + threshold;
    Bin_sizes(binInds) = Bins(b);
    Bin_offset(binInds) = offset-1;
    offset = offset + Bins(b);
end
offsets = ones(1,length(Bins));
for i = 2:length(Bins)
    offsets(i)  =  offsets(i) + sum(Bins(1:i-1));
end

%% Loop
est_spectrum = zeros(1,n);
wn = exp(1i*2*pi/n);
var_counter = 0;

for it = 1:k
    %% Take out zerotons
    signal_bins = (sum(abs(binSignals).^2,1)> T);
    Y = binSignals(:,signal_bins);
    n_signal_bins = sum(signal_bins);
    signal_2_actual_mapping = find(signal_bins) - Bin_offset(signal_bins);
    Bin_sizes_est = Bin_sizes(signal_bins);
    T_est         = T(signal_bins);
    
    %% Find singletons
    if (n_signal_bins)
        
        if (kay_flag)
            %% Estimate location
            est_location = zeros(1,n_signal_bins);
            for i = 1:B
                %%%% Kay starts here %%%
                M1 = Y( (i-1)*D + (1:D-1) , : );
                M2 = Y( (i-1)*D + (2:D) , : );
                
                omega = angle(conj(M1).*M2);
                % Find the omegas that are on boundary
                avg_high_freq = (median(abs(omega),1) > pi/2);
                omega  = omega + 2*pi*((omega < 0).*repmat(avg_high_freq,D-1,1));
                
                % Estimate locations using Kay's method
                wt_omega = wt*omega;
                wt_omega = wt_omega + 2*pi*(wt_omega < 0);
                est_location_b  =  mod(wt_omega*n/(2*pi),n);
                %%%% Kay ends here %%%
                
                % Update location (Take high order bit)
                nwrap = 2^(i-1);
                loc_update = est_location_b / nwrap - est_location;
                r =  loc_update - round(loc_update / (n/nwrap)) * (n/nwrap);
                est_location = est_location + r;
                
%               left = est_location - n/2^i;
%               right = est_location + n/2^i;
            end
            
            %% Process location estimates and estimate amplitude
            % Correcting the estimates using the windowing
            est_location = mod(round((est_location - (signal_2_actual_mapping - 1))./Bin_sizes_est).*Bin_sizes_est + (signal_2_actual_mapping - 1),n);
            
            % Estimate amplitude
            direction_matrix =  (wn.^(mod(delays.'*est_location,n)));
            est_alpha = sum(  Y.*conj(direction_matrix),1  )/N;
            est_energy = abs(est_alpha).^2;
            
            % Estimate the complete signal
            Signal_matrix = repmat(est_alpha,N,1).*direction_matrix;
            
            % Estimate noise matrix
            Noise_matrix  =   Y - Signal_matrix;
            noise_var_est = sum(abs(Noise_matrix).^2);
            % Set some of the singleton estimates to zero if noise is high
            singleton_est = (  noise_var_est <= T_est );
            singleton_est = singleton_est.*(est_energy > thresh);
            
            % Store to long arrays
            tmp = singleton_est;
            singleton_est = zeros(1,sum(Bins));
            singleton_est(signal_bins) = tmp;
            tmp = est_alpha;
            est_alpha = zeros(1,sum(Bins));
            est_alpha(signal_bins) = tmp;
            tmp = est_location;
            est_location = zeros(1,sum(Bins));
            est_location(signal_bins) = tmp;
            
        else
            singleton_est = zeros(1,sum(Bins));
            est_alpha = zeros(1,sum(Bins));
            est_location = zeros(1,sum(Bins));
            
            for stage = 1:length(Bins)
                n_bins = Bins(stage);
                for bin = 1:n_bins
                    minResid = inf;
                    for fold = 0:n/n_bins-1
                        bin_index = offsets(stage) + bin - 1;
                        l = (bin-1)+fold*n_bins;
                        freq = wn.^( delays'  * l);
                        alpha = mean(binSignals(:,bin_index) .* conj(freq));
                        sig = alpha*freq;
                        resid = sum(abs(binSignals(:,bin_index)-sig).^2,1);
                        
                        if (resid<minResid)
                            minResid = resid;
                            min_l = l;
                            min_alpha = alpha;
                        end
                    end
                    
                    if (minResid <= T(bin_index) && abs(min_alpha)^2 > thresh)
                        est_location(bin_index) = min_l;
                        est_alpha(bin_index) = min_alpha;
                        singleton_est(bin_index) = 1;
                    end
                end
            end
        end
        %         est_location
        %         est_location = est_location / B;
        
        %% Peel off singletons
        % Find all singletons and remove them
        singleton_found = 0;
        for stage = 1:length(Bins)
            n_bins = Bins(stage);
            for bin = 1:1:n_bins
                bin_index = offsets(stage) + bin - 1;
                is_singleton = singleton_est(bin_index);
                if(is_singleton)
                    est_val = est_alpha(bin_index);
                    est_loc = est_location(bin_index);
                    
                    if(mod(est_loc,n_bins) ~= (bin - 1))
                        continue;
                    end
                    
                    % Check if this is new singleton
                    if(  abs(est_spectrum(est_loc+1))==0   )
                        % At least one new singleton found
                        singleton_found = 1;
                        if(var_counter <= k)
                            est_spectrum(est_loc+1) = est_val;
                        end
                        var_counter = var_counter + 1;
                        
                        loc_signal = est_val.*(wn.^(mod(delays.'*est_loc,n)));
                        % Find the bin locations where it hashes
                        hash = mod(est_loc,Bins)+offsets;
                        binSignals(:,hash) = binSignals(:,hash) - repmat(loc_signal,1,length(hash));
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
est_xLoc = find(est_spectrum~=0);
