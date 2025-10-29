clear;
close all;
load sampleEEGdata

chan2use = 'fcz';
min_freq = 3;
max_freq = 30;
num_frex = 20;

chanidx = find(strcmpi({EEG.chanlocs.labels}, chan2use));
frex = logspace(log10(min_freq), log10(max_freq), num_frex);

fs = EEG.srate;
wave_number = 6;
nTrials = EEG.trials;
nTime = EEG.pnts;
nFreq = length(frex);
% freq * time * trials
power_all = zeros(nFreq, nTime, nTrials);

for tr = 1:nTrials
    signal = squeeze(EEG.data(chanidx, :, tr));
    [cfx, cfreqs] = morletWaveletTransform(signal, fs, frex, wave_number);
    power_all(:, :, tr) = abs(cfx).^2;
end
% avg on trials (freq*time)
avg_power=mean(power_all,3);

figure;
imagesc(EEG.times, frex, avg_power);
axis xy;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title(['Time-Frequency Power at ' upper(chan2use)]);
colorbar;



time_s = dsearchn(EEG.times',-500);
time_e = dsearchn(EEG.times',1200);
eegpower = power_all(:, time_s:time_e, :);
tftimes = EEG.times(time_s : time_e);
nTimepoints = length(tftimes);


%% parameters
 voxel_pval = 0.01;      % voxel-wise (uncorrected) treshold
 cluster_pval = 0.05;    % clutser-level alpha
 n_permutes = 2000;
base_idx = [dsearchn(tftimes',-500), dsearchn(tftimes',-100)] % baseline range

% Do a trial-level baseline normalization and store the trial-averaged
% baseline-normalized power in realmean
eegpower = power_all(:, time_s:time_e, :);  % got size as 20 x 436 x nTrials
nTimepoints = size(eegpower, 2);

% Baseline normalization
baseline_norm = zeros(nFreq, nTimepoints, nTrials);
for tr = 1:nTrials
    base_mean = mean(eegpower(:, base_idx(1):base_idx(2), tr), 2);
    baseline_norm(:, :, tr) = eegpower(:, :, tr) - base_mean;
end

% Trial-average
realmean = mean(baseline_norm, 3);  % got size as 20 x 436



%% shuffle the data  i.e; destroy the time locking to the stimuli to obtain 1000 equivalents of realmean and store in permuted vals

%% Permutation using circshift
nPerm = 1000;
permuted_vals = zeros(nFreq, nTimepoints, nPerm);

for p = 1:nPerm
    permuted_signal = zeros(nFreq, nTimepoints, nTrials);
    for tr = 1:nTrials
        shift_amt = randi(nTimepoints);
        for f = 1:nFreq
            permuted_signal(f, :, tr) = circshift(baseline_norm(f, :, tr), [0 shift_amt]);
        end
    end
    permuted_vals(:, :, p) = mean(permuted_signal, 3);  
end

%% create a z score metric 
perm_mean = mean(permuted_vals, 3);  % got size as 20 x 436
perm_std  = std(permuted_vals, [], 3);  % got size as 20 x 436

zmap = (realmean - perm_mean) ./ perm_std;  % got size as 20 x 436


%% cal those bins where z-score treshold exceeds pval (say , 0.05) .store results appropriately in treshmean
threshmean = realmean;
voxel_thresh = norminv(1-voxel_pval); % one tailed tresh
threshmean(abs(zmap) < voxel_thresh) = 0;


%% plot preliminary maps for realmean,zmap,and threshmean using contourf
figure;


subplot(1,3,1)
contourf(tftimes, frex, realmean, 40, 'linecolor', 'none');
axis xy;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Real Mean Power');
colorbar;


subplot(1,3,2)
contourf(tftimes, frex, zmap, 40, 'linecolor', 'none');
axis xy;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Z-map');
colorbar;

subplot(1,3,3)
contourf(tftimes, frex, threshmean, 40, 'linecolor', 'none');
axis xy;
xlabel('Time (ms)');
ylabel('Frequency (Hz)');
title('Thresholded Power (p < 0.01)');
colorbar;


colormap('jet');


%% FUNCTION: Morlet Wavelet Transform
function [cfx, cfreqso] = morletWaveletTransform(x, fs, cfreqs, morletParam, dim, plotFlag)
    if exist('morletParam', 'var') && strcmp(morletParam, 'plot')
        plotFlag = 'plot';
    end
    if ~exist('morletParam', 'var') || strcmp(morletParam, 'plot')
        morletParam = 7;
    end
    if ~isvector(x) && (~exist('dim', 'var') || strcmp(morletParam, 'plot'))
        dim = 2;
    end

    if ~isvector(x)
        permOrder = [dim, 1:dim-1, dim+1:ndims(x)];
        x = permute(x, permOrder);
        sx = size(x);
        x = x(:,:);
    end

    dt = 1/fs;
    morletFourierFactor = 4*pi/(morletParam + sqrt(2 + morletParam^2));
    scales = 1./(morletFourierFactor * cfreqs);

    if ~isvector(x)
        cfx = zeros([length(cfreqs), size(x)]);
    end

    if ~isvector(x)
        for ichan = 1:size(x,2)
            cfstruct = cwtft({x(:,ichan), dt}, 'scales', scales, 'wavelet', {'morl', morletParam});
            cfx(:, :, ichan) = cfstruct.cfs;
        end
        plotVal = squeeze(mean(cfx, 3));
    else
        cfstruct = cwtft({x, dt}, 'scales', scales, 'wavelet', {'morl', morletParam});
        cfx = cfstruct.cfs;
        plotVal = cfx;
    end

    sc = cfstruct.scales;
    cfreqso = 1./(sc * morletFourierFactor);

    if exist('plotFlag', 'var') && (strcmp(plotFlag, 'plot'))
        figure
        time = (1:length(x))/fs;
        imagesc(time, cfreqso, abs(plotVal))
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        title('Morlet wavelet power')
        axis xy

        figure
        imagesc(time, cfreqso, angle(plotVal))
        if exist('pmkmp_new', 'file') == 2
            colormap(pmkmp_new('ostwald_o', 256));
        else
            colormap(hsv)
        end
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        title('Morlet wavelet phase')
        axis xy
    end

    if ~isvector(x)
        cfx = reshape(cfx, [length(cfreqs), sx]);
        cfx = ipermute(cfx, [1, 1+permOrder]);
    end
end
