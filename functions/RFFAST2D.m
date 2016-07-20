function [est_spectrum,num_it] = sFFT2D_xy_random(k,noise_flag,kay_flag,binSignals,delay,BinsX,BinsY,nX,nY,thresh)

N = size(delay,1);
BinX_offset = zeros(sum(BinsX),max(BinsY));
BinX_sizes  = zeros(sum(BinsX),max(BinsY));
BinY_sizes  = zeros(sum(BinsX),max(BinsY));
T = zeros(sum(BinsX),max(BinsY));
threshold = 10^(-10);
noise_var   = (k./BinsX./BinsY)*noise_flag;
offsetX = 1;
for stage = 1:length(BinsX)
    binXInds = offsetX:(offsetX+BinsX(stage)-1);
    binYInds = 1:BinsY(stage);
    
    T(binXInds,binYInds) = 4*N*noise_var(stage) + threshold;
    BinX_sizes(binXInds,binYInds) = BinsX(stage);
    BinY_sizes(binXInds,binYInds) = BinsY(stage);
    BinX_offset(binXInds,binYInds) = offsetX-1;
    
    offsetX = offsetX + BinsX(stage);
end
offsetsX = ones(1,length(BinsX));
for i = 2:1:length(BinsX)
    offsetsX(i)  =  offsetsX(i) + sum(BinsX(1:i-1));
end

%%
est_spectrum = zeros(nX,nY);
wnX = exp(1i*2*pi/nX);
wnY = exp(1i*2*pi/nY);


num_singleton = 0;
num_it = 0;



est_alpha = zeros(sum(BinsX)*max(BinsY),1);
est_location = zeros(sum(BinsX)*max(BinsY),2);
est_stage = zeros(sum(BinsX)*max(BinsY),2);

while(1)
    
    % Find Singletons
    signal_bins = (squeeze(sum(abs(binSignals).^2,1))> T);
    
    if (sum(sum(signal_bins))==0)
        return;
    end
    
    % Map: check singleton and estimate alpha
    singleton_ind = 1;
    for stage = 1:1:length(BinsX)
        n_bins = BinsX(stage);
        for binX = 1:1:n_bins
            for binY = 1:1:n_bins
                bx = offsetsX(stage) + binX - 1;
                by = binY;
                
                minResid = inf;
                if (signal_bins(bx,by))
                    if (~kay_flag)
                        for foldX = 0:nX/n_bins-1
                            for foldY = 0:nY/n_bins-1
                                lx = (binX-1)+foldX*n_bins;
                                ly = (binY-1)+foldY*n_bins;
                                freq = wnX.^( delay(:,1,stage)  *lx).*wnY.^( delay(:,2,stage) *ly);
                                alpha = mean(binSignals(:,bx,by).*conj(freq));
                                sig = alpha*freq;
                                resid = sum(abs(binSignals(:,bx,by)-sig).^2,1);
                                
                                if (resid<minResid)
                                    minResid = resid;
                                    min_lx = lx;
                                    min_ly = ly;
                                    min_alpha = alpha;
                                end
                            end
                        end
                        
                        if (minResid <= T(bx,by) && abs(min_alpha)^2 > thresh)
                            est_location(singleton_ind,:) = [min_lx,min_ly];
                            est_alpha(singleton_ind) = min_alpha;
                            est_stage(singleton_ind) = stage;
                            singleton_ind = singleton_ind + 1;
                        end
                    else
                        % Kay's Single Frequency Estimate
                        % #### X ####
                        M1x = binSignals(1:N-1,bx,by);
                        M2x = binSignals(2:N,bx,by);
                        
                        omegaX = angle(conj(M1x).*M2x);
                        % find the omegas that are on boundary
                        avg_positive = (median(abs(omegaX),1) > pi/2);
                        omegaX  = omegaX + 2*pi*((omegaX < 0).*repmat(avg_positive,N-1,1));
                        
                        % apply weights
                        tX  = 0:N-2;
                        wtX = (1.5*N/(N^2-1))*( 1 - (  (tX-(N/2-1))/(N/2)  ).^2   );
                        
                        wt_omegaX = wtX*omegaX;
                        wt_omegaX = wt_omegaX + 2*pi*(wt_omegaX < 0);
                        lx  =  mod(round(wt_omegaX*nX/(2*pi)),nX);
                        % correcting the estimates using the windowing
                        lx = mod(round((lx - (bx - BinX_offset(bx,by) - 1))./BinX_sizes(bx,by)).*BinX_sizes(bx,by) + (bx - BinX_offset(bx,by) - 1),nX);
                        
                        freq = wnX.^( mod (delay(:,stage)  *lx , nX) );
                        alpha = mean(binSignals(:,bx,by).*conj(freq));
                        sig = alpha*freq;
                        resid = sum(abs(binSignals(:,bx,by)-sig).^2,1);
                        
                        % #### Y ####
                        M1y = binSignals(1:N-1,bx,by);
                        M2y = binSignals(2:N,bx,by);
                        
                        omegaY = angle(conj(M1y).*M2y);
                        % find the omegas that are on boundary
                        avg_positive = (median(abs(omegaY),1) > pi/2);
                        omegaY  = omegaY + 2*pi*((omegaY < 0).*repmat(avg_positive,N-1,1));
                        
                        % apply weights
                        tY  = 0:N-2;
                        wtY = (1.5*N/(N^2-1))*( 1 - (  (tY-(N/2-1))/(N/2)  ).^2   );
                        
                        wt_omegaY = wtY*omegaY;
                        wt_omegaY = wt_omegaY + 2*pi*(wt_omegaY < 0);
                        ly  =  mod(round(wt_omegaY*nY/(2*pi)),nY);
                        % correcting the estimates using the windowing
                        ly = mod(round((ly - (by - 1))./BinY_sizes(bx,by)).*BinY_sizes(bx,by) + (by - 1),nY);
                        
                        freqY = wnY.^( mod ( delay(:,stage) *ly , nY) );
                        alphaY = mean(binSignals(:,bx,by).*conj(freqY));
                        sigY = alphaY*freqY;
                        residY = sum(abs(binSignals(:,bx,by)-sigY).^2,1);
                        resid = resid+residY;
                        alpha = (alpha+alphaY)/2;
                        
                        if ((resid <= T(bx,by)) && (abs(alpha)^2 > thresh))
                            est_location(singleton_ind,:) = [lx,ly];
                            est_alpha(singleton_ind) = alpha;
                            est_stage(singleton_ind) = stage;
                            singleton_ind = singleton_ind + 1;
                        end
                    end
                    
                end
            end
        end
    end
    
    num_singleton_it = singleton_ind-1;
    
    
    if num_singleton_it==0
        return;
    end
    
    % Reduce: Put to spectrum, throw away invalid singletons
    num_singleton_it2 = 0;
%     est_spectrum_it2 = zeros(nX,nY);
    for singleton_ind = 1:num_singleton_it
        lx = est_location(singleton_ind,1);
        ly = est_location(singleton_ind,2);
        if (est_spectrum(lx+1,ly+1)==0)
            alpha = est_alpha(singleton_ind);
            
            
            invalidSingleton = 0;
            if (noise_flag)
                for stage = 1:length(BinsX)
                    % Find the bin locations where it hashes
                    hashX = mod(lx,BinsX(stage))+offsetsX(stage);
                    hashY = mod(ly,BinsY(stage))+1;
                    if ( stage ~= est_stage(singleton_ind) && (sum(binSignals(:,hashX,hashY).^2)==0))
                        invalidSingleton = 1;
                        break;
                    end
                end
            end
            
            if (~invalidSingleton)
                
                for stage = 1:length(BinsX)
                    
                    sig = alpha.*(wnX.^( mod ( delay(:,1,stage) , nX / BinsX(stage) )*lx)).*(wnY.^( mod ( delay(:,2,stage) , nY / BinsY(stage) )*ly));
                    
                    % Find the bin locations where it hashes
                    hashX = mod(lx,BinsX(stage))+offsetsX(stage);
                    hashY = mod(ly,BinsY(stage))+1;
                    binSignals(:,hashX,hashY) = binSignals(:,hashX,hashY) - sig;
                    if ( stage == est_stage(singleton_ind) )
                        binSignals(:,hashX,hashY) = 0;
                        binSignals(:,hashX,hashY) = 0;
                    end
                end
                
                est_spectrum(lx+1,ly+1) = alpha;
                
                num_singleton_it2 = num_singleton_it2 + 1;
            end
            
        end
    end
    
    if num_singleton_it2 == 0
        return;
    end
    
    num_it = num_it + 1;
    
    num_singleton = num_singleton + num_singleton_it2;
    if num_singleton > k
        return;
    end
    
    
end


