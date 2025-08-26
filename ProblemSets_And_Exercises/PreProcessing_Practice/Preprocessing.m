% -------------------------
% Load EEG dataset
% -------------------------
eeglab; 
EEG = pop_loadset('eeglab_data.set', fullfile(fileparts(which('eeglab')), 'sample_data'));

data = EEG.data;      
chan_labels = {EEG.chanlocs.labels}; 

% -------------------------
% Step 1: Filtering: Notch filter at 60 Hz
% -------------------------
wo = 60/(EEG.srate/2);   % Normalized notch frequency (w0 = 60Hz / Nyquist frequency)
bo = wo/35;              % Bandwidth parameter (controls how wide the notch is)  


% - 'wo' sets the center frequency of the notch filter.
% - 'bo' sets the bandwidth of the notch (how sharp or wide the stop-band is).
% If we replace "35" with a smaller number => filter becomes *wider* (removes more frequencies but may distort signal).
% If we replace "35" with a larger number => filter becomes *narrower* (more selective, less distortion but might leave some 60 Hz noise.).


[bn,an] = iirnotch(wo, bo);  

% Zero-phase filtering using filtfilt:
% - Normal Filtering Shifts Phase(Signal Delay in time)
% - Applies filter forward & backward to avoid phase distortion
% - keeps EEG peaks aligned in time (important for ERP analysis).
data_filt = filtfilt(bn, an, data')'; 

% -------------------------
% Step 2: Plot PSD (Power Spectral Density) before & after filtering
% -------------------------

spectopo(data, size(data,2), EEG.srate); % Before filtering
% Useful to visualize line noise (spikes) and check filter effects.
figure;
spectopo(data_filt, size(data_filt,2), EEG.srate); % After filtering

% -------------------------
% Step 3: Referencing: Common Average Reference
% -------------------------
% Referencing removes common noise across electrodes
% Each channel signal is re-referenced by subtracting the average across all channels
avg_ref = mean(data_filt, 1);
data_car = data_filt - avg_ref;

% -------------------------
% Step 4: Epoching around "rt" events
% -------------------------

event_latencies = [EEG.event.latency]; 
event_types     = {EEG.event.type};
target_idx = find(strcmp(event_types,'rt')); % Find "rt" events

epoch_window = round([-0.2 0.8]*EEG.srate);  % -200ms to +800ms window
epoch_len = diff(epoch_window)+1;

% Preallocate epochs: [channels x samples per epoch x number of events]
epochs = nan(size(data_car,1), epoch_len, length(target_idx));

for i = 1:length(target_idx)
    center = round(event_latencies(target_idx(i))); % Event latency in samples
    idx = center+epoch_window(1):center+epoch_window(2);
    if idx(1)>0 && idx(end)<=size(data_car,2)
        epochs(:,:,i) = data_car(:,idx);
    end
end

% -------------------------
% Step 5: Baseline correction
% -------------------------
% Baseline correction removes slow drifts by subtracting pre-stimulus mean activity
baseline_idx = 1:round(0.2*EEG.srate); % First 200ms = baseline window
baseline = mean(epochs(:,baseline_idx,:),2);
epochs_bc = epochs - baseline; % Trial-level baseline correction

% Global baseline correction: instead of trial-by-trial, subtract one global baseline
% Global baseline ensures consistency across trials but may mask trial variability.
global_baseline = mean(epochs(:,baseline_idx,:), [2 3]); % mean across channels & trials
epochs_gbc = epochs - global_baseline;

% -------------------------
% Step 5b: Time-Frequency Analysis
% -------------------------
% Compare trial- and channel-averaged spectrograms for epochs_bc vs epochs_gbc
fs    = EEG.srate;                    % sampling rate
win   = hamming(round(0.5*fs));       % 500 ms window
nover = round(0.4*fs);                % 400 ms overlap
nfft  = 2^nextpow2(length(win));      % FFT length

nChans  = size(epochs_bc,1);
nTrials = size(epochs_bc,3);

% === Epochs with trial-level baseline correction ===
P_all = [];
for ch = 1:nChans
    for tr = 1:nTrials
        sig = squeeze(epochs_bc(ch,:,tr));
        [~,F,T,P] = spectrogram(sig, win, nover, nfft, fs);
        if isempty(P_all), P_all = zeros([size(P), nChans, nTrials]); end
        P_all(:,:,ch,tr) = P;
    end
end
P_mean_bc = mean(P_all, [3 4]);

% === Epochs with global baseline correction ===
P_all = [];
for ch = 1:nChans
    for tr = 1:nTrials
        sig = squeeze(epochs_gbc(ch,:,tr));
        [~,F,T,P] = spectrogram(sig, win, nover, nfft, fs);
        if isempty(P_all), P_all = zeros([size(P), nChans, nTrials]); end
        P_all(:,:,ch,tr) = P;
    end
end
P_mean_gbc = mean(P_all, [3 4]);

% === Plot Trial- and Channel-averaged Spectrograms ===
figure;
surf(T, F, 10*log10(P_mean_bc), 'EdgeColor', 'none');
axis tight; view(0,90);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Trial- & Channel-averaged Spectrogram (Baseline-corrected)');
colorbar;

figure;
surf(T, F, 10*log10(P_mean_gbc), 'EdgeColor', 'none');
axis tight; view(0,90);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Trial- & Channel-averaged Spectrogram (Global Baseline)');
colorbar;

%% -------------------------
% Step 6: Improved Artifact Rejection
% -------------------------
% EEG artifacts = eye blinks, muscle noise, bad channels, etc.
% Rejecting trials and channels with abnormally high variance

% (a) Reject bad trials
trial_var = squeeze(var(epochs_bc,0,2)); % variance per channel & trial
trial_reject = (mean(trial_var,1) > mean(mean(trial_var))+3*std(mean(trial_var)));
good_trials = ~trial_reject;
epochs_clean = epochs_bc(:,:,good_trials);
fprintf('Rejected %d/%d trials\n', sum(trial_reject), length(trial_reject));

% (b) Reject bad channels
chan_var = squeeze(var(epochs_clean,0,[2 3])); % variance per channel
bad_chans = (chan_var < 1e-6) | (chan_var > mean(chan_var)+3*std(chan_var));
good_chans = find(~bad_chans);
epochs_clean = epochs_clean(good_chans,:,:);
fprintf('Rejected %d/%d channels\n', sum(bad_chans), length(chan_labels));

chan_labels_clean = chan_labels(good_chans);
