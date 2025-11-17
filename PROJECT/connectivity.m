%% Initialize EEGLAB
eeglab;
close all;

%% Define event mappings (equivalent to Python event_id dictionary)
event_id = containers.Map();
event_id('startOfNotRecognisedClip') = 1;
event_id('startOfRecognisedClipFirstWatch') = 2;
event_id('startOfRememberedClipFirstWatch') = 3;
event_id('startOfSecondWatch') = 128;
event_id('endOfClip') = 129;
event_id('recognitionClick') = 160;
event_id('trialStart') = 192;
event_id('trialEnd') = 193;
event_id('firstRecognitionClick') = 257;

%% Dataset parameters
participant = 1; 
bids_root = '.'; 
subject = sprintf('sub-%03d', participant);
task = 'MovieMemory';
suffix = 'eeg';
datatype = 'eeg';

%% Construct BIDS path
bids_path = fullfile(bids_root, subject, datatype, ...
    sprintf('%s_task-%s_%s.edf', subject, task, suffix));

fprintf('=== LOADING RAW DATA ===\n');
fprintf('Loading data from: %s\n', bids_path);

%% Load events first
events_file = fullfile(bids_root, subject, datatype, ...
    sprintf('%s_task-%s_events.tsv', subject, task));

fprintf('Loading events from: %s\n', events_file);
events_table = readtable(events_file, 'FileType', 'text', 'Delimiter', '\t');
fprintf('Found %d events in BIDS file\n', height(events_table));


%% Load EEG data

EEG = pop_biosig(bids_path, 'channels', 1:64);
EEG.setname = sprintf('sub-%03d_MovieMemory_raw', participant);
EEG = eeg_checkset(EEG);

fprintf('Successfully loaded EEG data:\n');
fprintf('  Subject: %s\n', subject);
fprintf('  Sampling rate: %d Hz\n', EEG.srate);
fprintf('  Number of channels: %d\n', EEG.nbchan);
fprintf('  Duration: %.2f seconds (%.2f minutes)\n', EEG.pnts/EEG.srate, (EEG.pnts/EEG.srate)/60);
fprintf('  Data points per channel: %d\n', EEG.pnts);


%% Replace events with BIDS events
EEG.event = [];
for i = 1:height(events_table)
    EEG.event(i).latency = events_table.onset(i) * EEG.srate + 1;
    EEG.event(i).type = events_table.trial_type{i};
    if ismember('duration', events_table.Properties.VariableNames)
        EEG.event(i).duration = events_table.duration(i) * EEG.srate;
    end
    if isKey(event_id, events_table.trial_type{i})
        EEG.event(i).code = event_id(events_table.trial_type{i});
    else
        EEG.event(i).code = 0;
    end
end
fprintf('Loaded %d events from BIDS file\n', length(EEG.event));

%% Visualization 1: Raw Data Overview
fprintf('\n=== VISUALIZING RAW DATA ===\n');
fprintf('Opening raw data visualization...\n'); % Raw Data - can have large spikes, flat lines, or excessive noise.

figure('Name', 'Raw Data - First 30 seconds', 'Position', [100 100 1200 800]);
plot_timerange = [0 30]; % First 30 seconds
time_samples = plot_timerange * EEG.srate + 1;
plot_data = EEG.data(:, time_samples(1):time_samples(2));
time_axis = (0:size(plot_data,2)-1) / EEG.srate;

channels_to_plot = 1:8:64; % Every 8th channel
for i = 1:length(channels_to_plot)
    ch = channels_to_plot(i);
    subplot(length(channels_to_plot), 1, i);
    plot(time_axis, plot_data(ch, :) + i*100); % Offset for visibility
    ylabel(sprintf('Ch%d', ch));
    if i == 1
        title('Raw EEG Data - First 30 seconds (Every 8th Channel)');
    end
    if i == length(channels_to_plot)
        xlabel('Time (seconds)');
    end
end

pop_eegplot(EEG, 1, 1, 1);

%% Calculate and display data statistics
fprintf('\n=== RAW DATA STATISTICS ===\n');
data_std = std(EEG.data, 0, 2);
data_mean = mean(EEG.data, 2);
fprintf('Channel standard deviations: Mean=%.2f μV, Range=[%.2f, %.2f] μV\n', ...
    mean(data_std), min(data_std), max(data_std));

% Check for flat channels
flat_channels = find(data_std < 1); % Less than 1 μV std
if ~isempty(flat_channels)
    fprintf('Potentially flat channels: %s\n', ...
        strjoin(string(flat_channels), ', '));
end

%% Step 1: Downsampling
fprintf('\n=== STEP 1: DOWNSAMPLING ===\n'); % From 2048Hz to 512Hz - Easy to work with computationally

original_srate = EEG.srate;
EEG = pop_resample(EEG, 512);
EEG.setname = 'Downsampled_512Hz';
EEG = eeg_checkset(EEG);

fprintf('Data points reduced from %d to %d per channel\n', ...
    EEG.pnts * (original_srate/512), EEG.pnts);

%% Visualization 2: Effect of Downsampling
figure('Name', 'Effect of Downsampling', 'Position', [200 200 1200 600]);
sample_channel = 32; % Fz approximately
time_segment = 1:min(5*original_srate, size(EEG.data,2)); % First 5 seconds

subplot(1,1,1);
original_time = (0:length(time_segment)-1) / original_srate;

plot((0:EEG.pnts-1)/EEG.srate, EEG.data(sample_channel,:));
title(sprintf('Channel %d After Downsampling (512 Hz)', sample_channel));
xlabel('Time (s)'); ylabel('Amplitude (μV)');


%% Step 2: High Pass Filtering
fprintf('\n=== STEP 2: HIGH PASS FILTERING ===\n'); % To remove drifts, Cutoff - 0.1Hz

data_before_hp = EEG.data;

EEG = pop_eegfiltnew(EEG, 0.1, [], [], false, [], 0);
EEG.setname = 'High_Pass_Filtered';
EEG = eeg_checkset(EEG);

fprintf('High-pass filtering complete.\n');

%% Visualization 3: Effect of High-Pass Filtering
figure('Name', 'Effect of High-Pass Filtering', 'Position', [300 300 1200 800]);

sample_channel = 32;
time_axis = (0:EEG.pnts-1) / EEG.srate;

subplot(3,1,1);
plot(time_axis, data_before_hp(sample_channel,:));
title(sprintf('Channel %d Before High-Pass Filter', sample_channel));
ylabel('Amplitude (μV)');

subplot(3,1,2);
plot(time_axis, EEG.data(sample_channel,:));
title(sprintf('Channel %d After High-Pass Filter (0.1 Hz)', sample_channel));
ylabel('Amplitude (μV)');

subplot(3,1,3);
plot(time_axis, data_before_hp(sample_channel,:) - EEG.data(sample_channel,:));
title('Removed Components (Low Frequency Drift)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');

%% Step 3: Bad Channel Detection and Removal
fprintf('\n=== STEP 3: BAD CHANNEL DETECTION ===\n'); % Using pop_clean_rawdata, channels flat for >5s, channels correlated <0.8 with neighbors, excessive 50/60 Hz noise 

originalEEG = EEG;
original_channels = {originalEEG.chanlocs.labels};

EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion', 5, 'ChannelCriterion', 0.8, ...
                       'LineNoiseCriterion', 4, 'Highpass', 'off', ...
                       'BurstCriterion', 'off', 'WindowCriterion', 'off');

remaining_channels = {EEG.chanlocs.labels};
removed_channels = setdiff(original_channels, remaining_channels);

fprintf('Bad channel detection complete:\n');
fprintf('  Original channels: %d\n', length(original_channels));
fprintf('  Remaining channels: %d\n', length(remaining_channels));
fprintf('  Removed channels: %s\n', strjoin(removed_channels, ', '));

%% Visualization 4: Channel Quality Assessment
if ~isempty(removed_channels)
    figure('Name', 'Removed Bad Channels', 'Position', [400 400 1200 600]);
    
    removed_indices = [];
    for i = 1:length(removed_channels)
        idx = find(strcmp(original_channels, removed_channels{i}));
        if ~isempty(idx)
            removed_indices(end+1) = idx;
        end
    end
    
    if ~isempty(removed_indices)
        time_axis = (0:originalEEG.pnts-1) / originalEEG.srate;
        
        % Plot first few removed channels
        n_to_plot = min(4, length(removed_indices));
        for i = 1:n_to_plot
            subplot(n_to_plot, 1, i);
            ch_idx = removed_indices(i);
            plot(time_axis, originalEEG.data(ch_idx, :));
            title(sprintf('Removed Channel: %s', original_channels{ch_idx}));
            ylabel('Amplitude (μV)');
            if i == n_to_plot
                xlabel('Time (s)');
            end
        end
        sgtitle('Examples of Removed Bad Channels');
    end
end


fprintf('Opening interactive plot for manual inspection...\n');
pop_eegplot(EEG, 1, 1, 1);

%% Step 4: Average Referencing
fprintf('\n=== STEP 4: AVERAGE REFERENCING ===\n'); % Use average of all channels as reference
data_before_reref = EEG.data;

EEG = pop_reref(EEG, []);
EEG.setname = 'Average_Referenced';
EEG = eeg_checkset(EEG);

fprintf('Average referencing complete.\n');
fprintf('New reference: Average of %d channels\n', EEG.nbchan);

%% Visualization 5: Effect of Average Referencing
figure('Name', 'Effect of Average Referencing', 'Position', [500 500 1200 600]);

sample_channels = [8, 16, 32]; % Different regions
time_segment = 1:min(10*EEG.srate, EEG.pnts); % First 10 seconds
time_axis = (0:length(time_segment)-1) / EEG.srate;

for i = 1:length(sample_channels)
    subplot(length(sample_channels), 2, 2*i-1);
    plot(time_axis, data_before_reref(sample_channels(i), time_segment));
    title(sprintf('Ch%d Before Re-referencing', sample_channels(i)));
    ylabel('Amplitude (μV)');
    
    subplot(length(sample_channels), 2, 2*i);
    plot(time_axis, EEG.data(sample_channels(i), time_segment));
    title(sprintf('Ch%d After Average Reference', sample_channels(i)));
    ylabel('Amplitude (μV)');
    
    if i == length(sample_channels)
        subplot(length(sample_channels), 2, 2*i-1);
        xlabel('Time (s)');
        subplot(length(sample_channels), 2, 2*i);
        xlabel('Time (s)');
    end
end

%% Step 5: Low Pass Filtering
fprintf('\n=== STEP 5: LOW PASS FILTERING ===\n'); % 40Hz to remove muscle artifacts
data_before_lp = EEG.data;

EEG = pop_eegfiltnew(EEG, [], 40, [], false, [], 0);
EEG.setname = 'Band_Pass_Filtered';
EEG = eeg_checkset(EEG);

fprintf('Low-pass filtering complete (0.1-40 Hz bandpass).\n');

%% Visualization 6: Frequency Domain Analysis
figure('Name', 'Frequency Domain Analysis', 'Position', [600 600 1200 800]);

sample_channel = 32;

% Before low-pass
[pxx_before, f] = pwelch(data_before_lp(sample_channel,:), [], [], [], EEG.srate);
% After low-pass
[pxx_after, f] = pwelch(EEG.data(sample_channel,:), [], [], [], EEG.srate);

subplot(2,2,1);
semilogy(f, pxx_before);
title('Power Spectrum Before Low-Pass Filter');
xlabel('Frequency (Hz)'); ylabel('Power (μV²/Hz)');
xlim([0 100]);

subplot(2,2,2);
semilogy(f, pxx_after);
title('Power Spectrum After 40 Hz Low-Pass Filter');
xlabel('Frequency (Hz)'); ylabel('Power (μV²/Hz)');
xlim([0 100]);

subplot(2,2,3);
time_axis = (0:min(5*EEG.srate, EEG.pnts)-1) / EEG.srate;
plot(time_axis, data_before_lp(sample_channel, 1:length(time_axis)));
title('Time Domain Before Low-Pass');
xlabel('Time (s)'); ylabel('Amplitude (μV)');

subplot(2,2,4);
plot(time_axis, EEG.data(sample_channel, 1:length(time_axis)));
title('Time Domain After Low-Pass');
xlabel('Time (s)'); ylabel('Amplitude (μV)');

%% Step 6: Epoching
fprintf('\n=== STEP 6: EPOCHING ===\n'); % -1 to 3 -> Baseline 1s, Task 3sec

%  - startOfRememberedClipFirstWatch: Successfully remembered clips
%  - startOfRecognisedClipFirstWatch: Recognized but uncertain clips
%  - startOfNotRecognisedClip: Unrecognized clips

epoch_events = {'startOfRememberedClipFirstWatch', ...
                'startOfRecognisedClipFirstWatch', ...
                'startOfNotRecognisedClip'};

available_events = unique({EEG.event.type});
fprintf('Available event types: %s\n', strjoin(available_events, ', '));

EEG = pop_epoch(EEG, epoch_events, [-1 10], 'newname', 'Epoched_Data', ...
               'epochinfo', 'yes');
EEG = eeg_checkset(EEG);

fprintf('Epoching complete:\n');
fprintf('  Epoch length: 4 seconds (-1 to +3 s)\n');
fprintf('  Total epochs: %d\n', EEG.trials);
fprintf('  Sampling points per epoch: %d\n', EEG.pnts);

%% Step 7: Baseline Correction
fprintf('\n=== STEP 7: BASELINE CORRECTION ===\n'); % 200 because if we including everything -> there may be some effects from previous trial ending

data_before_baseline = EEG.data;

EEG = pop_rmbase(EEG, [-200 0]);
EEG.setname = 'Baseline_Corrected';
EEG = eeg_checkset(EEG);

fprintf('Baseline correction complete using [-200, 0] ms baseline.\n');

%% Visualization 7: Effect of Baseline Correction
figure('Name', 'Effect of Baseline Correction', 'Position', [700 700 1200 600]);

sample_channel = 32;
sample_epoch = 1;
time_axis = EEG.times;

subplot(2,1,1);
plot(time_axis, data_before_baseline(sample_channel, :, sample_epoch));
hold on;
plot([-200 0], [0 0], 'r-', 'LineWidth', 3);
title('Before Baseline Correction');
ylabel('Amplitude (μV)');
legend('EEG signal', 'Baseline period');

subplot(2,1,2);
plot(time_axis, EEG.data(sample_channel, :, sample_epoch));
hold on;
plot([-200 0], [0 0], 'r-', 'LineWidth', 3);
title('After Baseline Correction');
xlabel('Time (ms)'); ylabel('Amplitude (μV)');

%% Step 8: Artifact Rejection
fprintf('\n=== STEP 8: ARTIFACT REJECTION ===\n'); % Automatic method using pop_autorej

epochs_before = EEG.trials;

EEG = pop_autorej(EEG, 'nogui', 'on', 'threshold', 100, 'startprob', 5, ...
                  'maxrej', 5, 'eegplot', 'off');

epochs_after = EEG.trials;
epochs_rejected = epochs_before - epochs_after;

fprintf('Artifact rejection complete:\n');
fprintf('  Epochs before: %d\n', epochs_before);
fprintf('  Epochs after: %d\n', epochs_after);
fprintf('  Epochs rejected: %d (%.1f%%)\n', epochs_rejected, (epochs_rejected/epochs_before)*100);


fprintf('Opening final epoched data for inspection...\n');
pop_eegplot(EEG, 1, 1, 1);

EEG.setname = 'Final_Clean_Data';
EEG = eeg_checkset(EEG);

%% Final Visualizations and Analysis
fprintf('\n=== FINAL DATA ANALYSIS ===\n');

% Trial count analysis
if EEG.trials > 0
    epoch_types = {EEG.epoch.eventtype};
    remembered_trials = sum(strcmp(epoch_types, 'startOfRememberedClipFirstWatch'));
    recognised_trials = sum(strcmp(epoch_types, 'startOfRecognisedClipFirstWatch'));
    not_recognised_trials = sum(strcmp(epoch_types, 'startOfNotRecognisedClip'));
    
    % Visualization 8: Trial Distribution
    figure('Name', 'Trial Distribution by Condition', 'Position', [800 800 800 600]);
    
    subplot(2,2,1);
    trial_counts = [remembered_trials, recognised_trials, not_recognised_trials];
    trial_labels = {'Remembered', 'Recognised', 'Not Recognised'};
    bar(trial_counts);
    set(gca, 'XTickLabel', trial_labels);
    title('Number of Trials per Condition');
    ylabel('Number of Trials');
    
    subplot(2,2,2);
    pie(trial_counts, trial_labels);
    title('Trial Distribution (%)');
    
    % Sample ERPs
    subplot(2,2,[3,4]);
    sample_channel = 32; % Approximate Fz location
    
    % Calculate average ERPs for each condition
    if remembered_trials > 0
        remembered_indices = strcmp(epoch_types, 'startOfRememberedClipFirstWatch');
        erp_remembered = mean(EEG.data(sample_channel, :, remembered_indices), 3);
        plot(EEG.times, erp_remembered, 'b-', 'LineWidth', 2);
        hold on;
    end
    
    if recognised_trials > 0
        recognised_indices = strcmp(epoch_types, 'startOfRecognisedClipFirstWatch');
        erp_recognised = mean(EEG.data(sample_channel, :, recognised_indices), 3);
        plot(EEG.times, erp_recognised, 'g-', 'LineWidth', 2);
    end
    
    if not_recognised_trials > 0
        not_recognised_indices = strcmp(epoch_types, 'startOfNotRecognisedClip');
        erp_not_recognised = mean(EEG.data(sample_channel, :, not_recognised_indices), 3);
        plot(EEG.times, erp_not_recognised, 'r-', 'LineWidth', 2);
    end
    
    xlabel('Time (ms)');
    ylabel('Amplitude (μV)');
    title(sprintf('Grand Average ERPs - Channel %d', sample_channel));
    legend('Remembered', 'Recognised', 'Not Recognised');
    grid on;
    
else
    fprintf('No epochs were created - check event types and timing\n');
end

%% Final Summary
fprintf('\n=== PREPROCESSING SUMMARY ===\n');
fprintf('Original sampling rate: 2048 Hz\n');
fprintf('Final sampling rate: %d Hz\n', EEG.srate);
fprintf('Original channels: 64\n');
fprintf('Channels after cleaning: %d\n', EEG.nbchan);
if ~isempty(removed_channels)
    fprintf('Removed channels: %s\n', strjoin(removed_channels, ', '));
else
    fprintf('No channels were removed\n');
end
fprintf('Total epochs: %d\n', EEG.trials);
fprintf('Epoch length: %.1f seconds\n', (EEG.pnts-1)/EEG.srate);

if exist('remembered_trials', 'var')
    fprintf('\nTrial counts by condition:\n');
    fprintf('  Remembered clips: %d trials\n', remembered_trials);
    fprintf('  Recognised clips: %d trials\n', recognised_trials);
    fprintf('  Not recognised clips: %d trials\n', not_recognised_trials);
    fprintf('  Total: %d trials\n', remembered_trials + recognised_trials + not_recognised_trials);
    
    % Calculate recognition rates
    total_presented = remembered_trials + recognised_trials + not_recognised_trials;
    recognition_rate = (remembered_trials + recognised_trials) / total_presented * 100;
    memory_rate = remembered_trials / total_presented * 100;
    
    fprintf('\nBehavioral Performance:\n');
    fprintf('  Recognition rate: %.1f%%\n', recognition_rate);
    fprintf('  Strong memory rate: %.1f%%\n', memory_rate);
end

%% Save processed data
save_filename = sprintf('%s_task-%s_preprocessed.set', subject, task);
pop_saveset(EEG, 'filename', save_filename, 'filepath', pwd);
fprintf('\nPreprocessed data saved as: %s\n', save_filename);

fprintf('\n=== PREPROCESSING COMPLETE ===\n');
fprintf('Data is now ready for ERP analysis, time-frequency analysis, or connectivity analysis.\n');

%%
%% === SPLIT INTO ERP AND CONNECTIVITY DATASETS ===
ERP_EEG = EEG;   % will keep average reference for ERP
CONN_EEG = EEG;  % will apply bipolar reference for connectivity

%%
%% === DEFINE BIPOLAR PAIRS FOR CONNECTIVITY ===

bipolar_pairs = {
    % === LEFT ANTERIOR-POSTERIOR CHAIN ===
    'Fp1','AF7';
    'AF7','AF3';
    'AF3','F1';
    'F1','F3';
    'F3','F5';
    'F5','F7';
    'F7','FT7';
    'FT7','TP7';
    'TP7','CP5';
    'CP5','P7';
    'P7','PO7';
    'PO7','O1';

    % === RIGHT ANTERIOR-POSTERIOR CHAIN (corrected for your cap) ===
    'Fp2','AF8';
    'AF8','AF4';
    'AF4','F2';
    'F2','F4';
    'F4','F6';
    'F6','FC6';
    'FC6','C6';
    'C6','TP8';    % TP8 exists in your list
    'TP8','CP6';
    'CP6','P8';
    'P8','PO8';
    'PO8','O2';

    % === MIDLINE CHAIN (frontal -> parietal) ===
    'Fpz','AFz';
    'AFz','Fz';
    'Fz','FCz';
    'FCz','Cz';
    'Cz','CPz';
    'CPz','Pz';
    'Pz','POz';
    'POz','Oz';

    % === PARIETAL LATERAL SHORT CHAINS (extra coverage) ===
    'P3','P5';
    'P5','P7';
    'P4','P6';
    'P6','P8';
    'CP1','P3';    % centro-parietal short pair
    'CP2','P4';
};
%%

function EEG = apply_bipolar(EEG, pairs)

nPairs = size(pairs,1);
new_data = zeros(nPairs, EEG.pnts, EEG.trials);
new_chanlocs = struct([]);

for i = 1:nPairs
    ch1 = find(strcmpi({EEG.chanlocs.labels}, pairs{i,1}));
    ch2 = find(strcmpi({EEG.chanlocs.labels}, pairs{i,2}));

    if isempty(ch1) || isempty(ch2)
        fprintf('Skipping %s-%s (channel missing)\n', pairs{i,1}, pairs{i,2});
        continue;
    end
    
    % bipolar difference
    new_data(i,:,:) = EEG.data(ch1,:,:) - EEG.data(ch2,:,:);

    % create a new channel label
    new_chanlocs(i).labels = sprintf('%s-%s', pairs{i,1}, pairs{i,2});
end

% Create a NEW EEG structure for connectivity dataset
CONN = EEG;

% Replace data
CONN.data = new_data;
CONN.nbchan = size(new_data,1);

% Replace channel information
CONN.chanlocs = new_chanlocs;  % Only labels are needed for connectivity
CONN.urchanlocs = new_chanlocs;
CONN.chaninfo = [];

% Remove ICA weights (no longer valid)
CONN.icawinv = [];
CONN.icaweights = [];
CONN.icasphere = [];

% Update EEG structure
EEG = eeg_checkset(CONN);

fprintf('Bipolar referencing complete: %d bipolar channels created.\n', EEG.nbchan);
end

%%
CONN_EEG = apply_bipolar(ERP_EEG, bipolar_pairs);

CONN_EEG = eeg_checkset(CONN_EEG);


%%
%% ========= CONNECTIVITY ANALYSIS BLOCK ============
fprintf('\n=== SIMPLIFIED CONNECTIVITY ANALYSIS ===\n');
%% ========= STEP 0: BALANCE TRIALS ==========
fprintf('\n=== BALANCING TRIALS ACROSS CONDITIONS ===\n');

% Determine minimum trials across conditions
nTrials = zeros(1, length(conditions_to_test));
for i = 1:length(conditions_to_test)
    cond = conditions_to_test{i};
    if ~isempty(condition_data.(cond))
        nTrials(i) = condition_data.(cond).trials;
    else
        nTrials(i) = 0;
    end
end

minTrials = min(nTrials(nTrials>0));
fprintf('Minimum trials across conditions: %d\n', minTrials);

% Subsample each condition to minTrials
for i = 1:length(conditions_to_test)
    cond = conditions_to_test{i};
    if ~isempty(condition_data.(cond))
        idx = randperm(condition_data.(cond).trials, minTrials);
        condition_data.(cond).data = condition_data.(cond).data(:,:,idx);
        condition_data.(cond).trials = minTrials;
        fprintf('%s: reduced to %d trials\n', cond, minTrials);
    end
end


%% ========= STEP 1: DEFINE CONDITIONS ================
% 2-column cell array: first column = event name, second = label
conditions = {
    'startOfRememberedClipFirstWatch', 'Remembered';
    'startOfRecognisedClipFirstWatch', 'Recognised';
    'startOfNotRecognisedClip', 'Forgotten'
};

fprintf('Preparing condition datasets...\n\n');

%% ========= STEP 2: Extract data for each condition =============
condition_data = struct();

for i = 1:size(conditions,1)
    event_name = conditions{i,1};
    cond_name  = conditions{i,2};

    try
       temp_EEG = pop_selectevent(CONN_EEG, 'type', event_name);

        % Subsample to minTrials
        if temp_EEG.trials > minTrials
            idx = randperm(temp_EEG.trials, minTrials);
            temp_EEG.data = temp_EEG.data(:,:,idx);
            temp_EEG.trials = minTrials;
        end
        
        condition_data.(cond_name) = temp_EEG;


    catch
        fprintf('%s: NO trials found\n', cond_name);
        condition_data.(cond_name) = [];
    end
end



%% ========= STEP 3: Simple PLV Connectivity (theta 4-8 Hz) ==========
fprintf('\n=== COMPUTING SIMPLE THETA PLV (4–8 Hz) ===\n');

theta_results = struct();
conditions_to_test = {'Remembered','Recognised','Forgotten'};

for i = 1:length(conditions_to_test)
    cond = conditions_to_test{i};

    if ~isempty(condition_data.(cond))
        EEGtmp = condition_data.(cond);

        avg_plv = simple_plv_connectivity(EEGtmp, [4 8]);

        theta_results.(cond).avg_plv = avg_plv;
        theta_results.(cond).ntrials = EEGtmp.trials;

        fprintf('%s: PLV = %.3f (n=%d)\n', cond, avg_plv, EEGtmp.trials);
    else
        fprintf('%s: No data available for PLV computation\n', cond);
        theta_results.(cond).avg_plv = NaN;
        theta_results.(cond).ntrials = 0;
    end
end

%% ========= STEP 4: Basic Hypothesis Checks =============
fprintf('\n=== CONNECTIVITY COMPARISONS ===\n');

if all(isfield(theta_results, {'Remembered','Forgotten'}))
    d = theta_results.Remembered.avg_plv - theta_results.Forgotten.avg_plv;
    fprintf('Remembered − Forgotten = %.3f\n', d);
end

if all(isfield(theta_results, {'Remembered','Recognised'}))
    d = theta_results.Remembered.avg_plv - theta_results.Recognised.avg_plv;
    fprintf('Remembered − Recognised = %.3f\n', d);
end

if all(isfield(theta_results, {'Recognised','Forgotten'}))
    d = theta_results.Recognised.avg_plv - theta_results.Forgotten.avg_plv;
    fprintf('Recognised − Forgotten = %.3f\n', d);
end

%% ========= STEP 5: Visualization =======================
fprintf('\nPlotting results...\n');

figure('Position',[50 50 900 400]);

% Bar plot: average PLV
subplot(1,2,1);
conds = fieldnames(theta_results);
vals = cellfun(@(c) theta_results.(c).avg_plv, conds);
bar(vals);
set(gca,'XTick',1:length(conds),'XTickLabel',conds);
ylabel('Theta PLV');
title('Average Theta Connectivity');
grid on;

% Differences between conditions
subplot(1,2,2);
pairs = {'Rem-Forg','Rem-Rec','Rec-Forg'};
diffs = [];

% Remembered − Forgotten
if ~isnan(theta_results.Remembered.avg_plv) && ~isnan(theta_results.Forgotten.avg_plv)
    diffs(end+1) = theta_results.Remembered.avg_plv - theta_results.Forgotten.avg_plv;
else
    diffs(end+1) = NaN;
end

% Remembered − Recognised
if ~isnan(theta_results.Remembered.avg_plv) && ~isnan(theta_results.Recognised.avg_plv)
    diffs(end+1) = theta_results.Remembered.avg_plv - theta_results.Recognised.avg_plv;
else
    diffs(end+1) = NaN;
end

% Recognised − Forgotten
if ~isnan(theta_results.Recognised.avg_plv) && ~isnan(theta_results.Forgotten.avg_plv)
    diffs(end+1) = theta_results.Recognised.avg_plv - theta_results.Forgotten.avg_plv;
else
    diffs(end+1) = NaN;
end

bar(diffs);
set(gca,'XTick',1:3,'XTickLabel',pairs);
ylabel('PLV Difference');
title('Connectivity Differences');
grid on;

%% ========= SUPPORTING FUNCTIONS ================
function avg_plv = simple_plv_connectivity(EEG, band)
    data = EEG.data; % channels × time × trials
    [nch, ~, ntrials] = size(data);
    srate = EEG.srate;

    % Bandpass filter
    filtered = bandpass_filter(data, srate, band);

    % Initialize phase array
    phase = zeros(size(filtered));

    % Hilbert transform per channel and trial
    for ch = 1:nch
        for tr = 1:ntrials
            phase(ch,:,tr) = angle(hilbert(squeeze(filtered(ch,:,tr))));
        end
    end

    % Compute PLV across all pairs
    plv_all = [];
    for ch1 = 1:nch
        for ch2 = ch1+1:nch
            pdiff = phase(ch1,:,:) - phase(ch2,:,:);
            plv = abs(mean(exp(1i * pdiff), 'all'));
            plv_all(end+1) = plv;
        end
    end

    avg_plv = mean(plv_all);
end

function out = bandpass_filter(data, srate, band)
    % 4th order Butterworth bandpass filter applied to 3D EEG data
    % data: channels x time x trials
    [b,a] = butter(4, band/(srate/2), 'bandpass');
    [nch, ntime, ntrials] = size(data);
    out = zeros(nch, ntime, ntrials);

    for ch = 1:nch
        for tr = 1:ntrials
            out(ch,:,tr) = filtfilt(b, a, double(squeeze(data(ch,:,tr))));
        end
    end
end

%%
%% ========= STEP 6: COMPUTE OVERALL PLV ACROSS ALL CONDITIONS ============
fprintf('\n=== COMPUTING OVERALL THETA PLV ACROSS ALL CONDITIONS ===\n');

all_trials = [];
for i = 1:length(conditions_to_test)
    cond = conditions_to_test{i};
    if ~isempty(condition_data.(cond))
        all_trials = cat(3, all_trials, condition_data.(cond).data);
    end
end

if ~isempty(all_trials)
    overall_theta_PLV = simple_plv_connectivity(struct('data', all_trials, 'srate', CONN_EEG.srate), [4 8]);
    fprintf('Overall PLV (all conditions combined): %.3f\n', overall_theta_PLV);
else
    overall_theta_PLV = NaN;
    fprintf('No trials available across conditions for overall PLV\n');
end

%% ========= STEP 7: PERMUTATION TEST BETWEEN CONDITIONS ============
fprintf('\n=== PERMUTATION TEST BETWEEN CONDITIONS ===\n');

perm_pairs = {'Remembered','Forgotten'; 'Remembered','Recognised'; 'Recognised','Forgotten'};
nPerm = 1000; 
perm_results = struct();

for i = 1:size(perm_pairs,1)
    % if mod(p, 500) == 0  % print every 500 iterations to avoid too much output
    
   
    c1 = perm_pairs{i,1};
    c2 = perm_pairs{i,2};

    data1 = condition_data.(c1).data;
    data2 = condition_data.(c2).data;

    if isempty(data1) || isempty(data2)
        fprintf('%s vs %s: missing data, skipping\n', c1, c2);
        perm_results.(sprintf('%s_vs_%s', c1, c2)) = NaN;
        continue
    end

    % Observed PLV difference
    obs_diff = theta_results.(c1).avg_plv - theta_results.(c2).avg_plv;

    combined = cat(3, data1, data2);
    n1 = size(data1,3);
    n2 = size(data2,3);

    perm_diff = zeros(1, nPerm);
    for p = 1:nPerm
        fprintf('Permutation iteration %d / %d\n', p, nPerm);
        idx = randperm(n1+n2);
        perm_group1 = combined(:,:,idx(1:n1));
        perm_group2 = combined(:,:,idx(n1+1:end));

        plv1 = simple_plv_connectivity(struct('data', perm_group1,'srate',CONN_EEG.srate), [4 8]);
        plv2 = simple_plv_connectivity(struct('data', perm_group2,'srate',CONN_EEG.srate), [4 8]);

        perm_diff(p) = plv1 - plv2;
    end

    % Two-sided p-value
    p_val = mean(abs(perm_diff) >= abs(obs_diff));
    fprintf('%s vs %s: observed diff = %.3f, p = %.4f\n', c1, c2, obs_diff, p_val);

    perm_results.(sprintf('%s_vs_%s', c1, c2)).obs_diff = obs_diff;
    perm_results.(sprintf('%s_vs_%s', c1, c2)).p_val = p_val;
    perm_results.(sprintf('%s_vs_%s', c1, c2)).perm_diff = perm_diff;
end

%% ========= STEP 8: OPTIONAL: PLOT PERMUTATION DISTRIBUTIONS ============
figure('Position',[100 100 1200 400]);
for i = 1:size(perm_pairs,1)
    pair_name = sprintf('%s_vs_%s', perm_pairs{i,1}, perm_pairs{i,2});
    if isfield(perm_results, pair_name) && ~isnan(perm_results.(pair_name))
        subplot(1,3,i);
        histogram(perm_results.(pair_name).perm_diff,50);
        hold on;
        ylims = ylim;
        line([perm_results.(pair_name).obs_diff perm_results.(pair_name).obs_diff], [0 ylims(2)], ...
            'Color','r','LineWidth',2);
        title([perm_pairs{i,1} ' - ' perm_pairs{i,2}]);
        xlabel('PLV difference (permutation)');
        ylabel('Frequency');
        grid on;
    end
end
