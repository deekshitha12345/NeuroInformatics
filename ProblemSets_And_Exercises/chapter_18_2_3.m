

clear; clc; close all;

%% --- LOAD SAMPLE DATA ---
load sampleEEGdata  


data      = EEG.data;       
times     = EEG.times;      
chanlocs  = EEG.chanlocs;   
srate     = EEG.srate;     
[chanN, pnts, trialN] = size(data);

fprintf('Data: %d channels, %d time points, %d trials\n', chanN, pnts, trialN);

%% --- TIME-FREQUENCY PARAMETERS ---
frex     = 2:2:30;       
num_frex = length(frex);
n_cycles = 6;            


tfpower = zeros(chanN, num_frex, pnts);

%% --- TIME-FREQUENCY DECOMPOSITION ---
fprintf('Computing time–frequency decomposition...\n');

for ch = 1:chanN
  
    signal = mean(squeeze(data(ch,:,:)),2);
    

    n_wavelet = length(signal)*2;
    fft_signal = fft(signal, n_wavelet);
    
    for fi = 1:num_frex
        
        s = n_cycles / (2*pi*frex(fi));
        t = -2:1/srate:2;
        wavelet = exp(2*1i*pi*frex(fi).*t) .* exp(-t.^2./(2*s^2));
        fft_wavelet = fft(wavelet, n_wavelet);
        fft_wavelet = fft_wavelet ./ max(fft_wavelet);
        
       
        conv_res = ifft(fft_wavelet .* fft_signal);
        conv_res = conv_res(1:length(signal)); 
        
        % power
        tfpower(ch,fi,:) = abs(conv_res).^2;
    end
end

fprintf('TF decomposition done! Size of tfpower = %s\n', mat2str(size(tfpower)));

%% --- BASELINE NORMALIZATION AND TOPO MAPS ---
freqs_to_plot   = [6 10 20];       
timepoints_ms   = [100 200 300 400 500];
baseline_window = [-300 -100];

[~,baseidx(1)] = min(abs(times - baseline_window(1)));
[~,baseidx(2)] = min(abs(times - baseline_window(2)));
[~,timeidx]    = arrayfun(@(x) min(abs(times - x)), timepoints_ms);

for fi = 1:length(freqs_to_plot)
    
    
    [~,freq_idx] = min(abs(frex - freqs_to_plot(fi)));
    
  
    baseline_power = mean(tfpower(:,freq_idx,baseidx(1):baseidx(2)),3);
    
   
    tfpower_dB = 10 * log10(bsxfun(@rdivide, tfpower(:,freq_idx,:), baseline_power));
    
  
    clim_noBL = [-3 3];
    clim_BL   = [-5 5];
    
    figure('Name',['Freq ' num2str(frex(freq_idx)) ' Hz'],'Color','w');
    
    for ti = 1:length(timeidx)
        % no baseline normalization
        subplot(2,length(timeidx),ti);
        topoplot(tfpower(:,freq_idx,timeidx(ti)), chanlocs, 'maplimits', clim_noBL);
        title([num2str(times(timeidx(ti))) ' ms']);
        if ti==1, ylabel('No baseline'); end
        colorbar;
        
        % with baseline normalization
        subplot(2,length(timeidx),ti+length(timeidx));
        topoplot(tfpower_dB(:,timeidx(ti)), chanlocs, 'maplimits', clim_BL);
        title([num2str(times(timeidx(ti))) ' ms']);
        if ti==1, ylabel('With baseline'); end
        colorbar;
    end
    
    sgtitle(['Topographical Power Maps @ ' num2str(frex(freq_idx)) ' Hz']);
end
