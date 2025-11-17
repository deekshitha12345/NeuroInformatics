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
fprintf('Opening raw data visualization...\n');
figure('Name', 'Raw Data - First 30 seconds', 'Position', [100 100 1200 800]);
plot_timerange = [0 30];
time_samples = plot_timerange * EEG.srate + 1;
plot_data = EEG.data(:, time_samples(1):time_samples(2));
time_axis = (0:size(plot_data,2)-1) / EEG.srate;
channels_to_plot = 1:8:64;
for i = 1:length(channels_to_plot)
    ch = channels_to_plot(i);
    subplot(length(channels_to_plot), 1, i);
    plot(time_axis, plot_data(ch, :) + i*100);
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
flat_channels = find(data_std < 1);
if ~isempty(flat_channels)
    fprintf('Potentially flat channels: %s\n', ...
        strjoin(string(flat_channels), ', '));
end

%% Step 1: Downsampling
fprintf('\n=== STEP 1: DOWNSAMPLING ===\n');
original_srate = EEG.srate;
EEG = pop_resample(EEG, 512);
EEG.setname = 'Downsampled_512Hz';
EEG = eeg_checkset(EEG);
fprintf('Data points reduced from %d to %d per channel\n', ...
    EEG.pnts * (original_srate/512), EEG.pnts);

%% Visualization 2: Effect of Downsampling
figure('Name', 'Effect of Downsampling', 'Position', [200 200 1200 600]);
sample_channel = 32;
time_segment = 1:min(5*original_srate, size(EEG.data,2));
subplot(1,1,1);
plot((0:EEG.pnts-1)/EEG.srate, EEG.data(sample_channel,:));
title(sprintf('Channel %d After Downsampling (512 Hz)', sample_channel));
xlabel('Time (s)'); ylabel('Amplitude (μV)');

%% Step 2: High Pass Filtering
fprintf('\n=== STEP 2: HIGH PASS FILTERING ===\n');
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
fprintf('\n=== STEP 3: BAD CHANNEL DETECTION ===\n');
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
fprintf('\n=== STEP 4: AVERAGE REFERENCING ===\n');
data_before_reref = EEG.data;
EEG = pop_reref(EEG, []);
EEG.setname = 'Average_Referenced';
EEG = eeg_checkset(EEG);
fprintf('Average referencing complete.\n');
fprintf('New reference: Average of %d channels\n', EEG.nbchan);

%% Visualization 5: Effect of Average Referencing
figure('Name', 'Effect of Average Referencing', 'Position', [500 500 1200 600]);
sample_channels = [8, 16, 32];
time_segment = 1:min(10*EEG.srate, EEG.pnts);
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
fprintf('\n=== STEP 5: LOW PASS FILTERING ===\n');
data_before_lp = EEG.data;
EEG = pop_eegfiltnew(EEG, [], 40, [], false, [], 0);
EEG.setname = 'Band_Pass_Filtered';
EEG = eeg_checkset(EEG);
fprintf('Low-pass filtering complete (0.1-40 Hz bandpass).\n');

%% Visualization 6: Frequency Domain Analysis
figure('Name', 'Frequency Domain Analysis', 'Position', [600 600 1200 800]);
sample_channel = 32;
[pxx_before, f] = pwelch(data_before_lp(sample_channel,:), [], [], [], EEG.srate);
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
fprintf('\n=== STEP 6: EPOCHING ===\n');
epoch_events = {'startOfRememberedClipFirstWatch', ...
                'startOfRecognisedClipFirstWatch', ...
                'startOfNotRecognisedClip'};
available_events = unique({EEG.event.type});
fprintf('Available event types: %s\n', strjoin(available_events, ', '));
EEG = pop_epoch(EEG, epoch_events, [-1 3], 'newname', 'Epoched_Data', ...
               'epochinfo', 'yes');
EEG = eeg_checkset(EEG);
fprintf('Epoching complete:\n');
fprintf('  Epoch length: 4 seconds (-1 to +3 s)\n');
fprintf('  Total epochs: %d\n', EEG.trials);
fprintf('  Sampling points per epoch: %d\n', EEG.pnts);

%% Step 7: Baseline Correction
fprintf('\n=== STEP 7: BASELINE CORRECTION ===\n');
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
fprintf('\n=== STEP 8: ARTIFACT REJECTION ===\n');
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
if EEG.trials > 0
    epoch_types = {EEG.epoch.eventtype};
    remembered_trials = sum(strcmp(epoch_types, 'startOfRememberedClipFirstWatch'));
    recognised_trials = sum(strcmp(epoch_types, 'startOfRecognisedClipFirstWatch'));
    not_recognised_trials = sum(strcmp(epoch_types, 'startOfNotRecognisedClip'));
    
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
    subplot(2,2,[3,4]);
    sample_channel = 32;
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
ERP_EEG = EEG;
CONN_EEG = EEG;

%%
%% === DEFINE BIPOLAR PAIRS FOR CONNECTIVITY (Full Montage) ===
bipolar_pairs = {
    % === LEFT ANTERIOR-POSTERIOR CHAIN ===
    'Fp1','AF7'; 'AF7','AF3'; 'AF3','F1'; 'F1','F3'; 'F3','F5'; 'F5','F7'; 'F7','FT7'; 
    'FT7','TP7'; 'TP7','CP5'; 'CP5','P7'; 'P7','PO7'; 'PO7','O1';
    % === RIGHT ANTERIOR-POSTERIOR CHAIN (corrected for your cap) ===
    'Fp2','AF8'; 'AF8','AF4'; 'AF4','F2'; 'F2','F4'; 'F4','F6'; 'F6','FC6'; 'FC6','C6'; 
    'C6','TP8'; 'TP8','CP6'; 'CP6','P8'; 'P8','PO8'; 'PO8','O2';
    % === MIDLINE CHAIN (frontal -> parietal) ===
    'Fpz','AFz'; 'AFz','Fz'; 'Fz','FCz'; 'FCz','Cz'; 'Cz','CPz'; 'CPz','Pz'; 'Pz','POz'; 'POz','Oz';
    % === PARIETAL LATERAL SHORT CHAINS (extra coverage) ===
    'P3','P5'; 'P5','P7'; 'P4','P6'; 'P6','P8'; 'CP1','P3'; 'CP2','P4';
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
    
    new_data(i,:,:) = EEG.data(ch1,:,:) - EEG.data(ch2,:,:);
    new_chanlocs(i).labels = sprintf('%s-%s', pairs{i,1}, pairs{i,2});
end

CONN = EEG;
CONN.data = new_data;
CONN.nbchan = size(new_data,1);
CONN.chanlocs = new_chanlocs;
CONN.urchanlocs = new_chanlocs;
CONN.chaninfo = [];
CONN.icawinv = []; CONN.icaweights = []; CONN.icasphere = [];
EEG = eeg_checkset(CONN);
fprintf('Bipolar referencing complete: %d bipolar channels created.\n', EEG.nbchan);
end

%% Execute Bipolar Referencing
CONN_EEG = apply_bipolar(ERP_EEG, bipolar_pairs);
CONN_EEG = eeg_checkset(CONN_EEG);

%%
%% === DEFINE TARGET BIPOLAR PAIRS (Hypothesized Network Groups) ===
function [frontal_nodes, posterior_nodes] = get_target_network_groups()
    % Frontal/Prefrontal Nodes (F) - Represents PFC/Frontal Lobe Activity
    frontal_nodes = {
        'Fp1-AF7', 'AF7-AF3', 'AF3-F1', 'F1-F3', 'F3-F5', 'F5-F7', 'Fz-FCz', 'AFz-Fz', 'Fp2-AF8', 'AF8-AF4', 'AF4-F2', 'F2-F4', 'F4-F6', 'F6-FC6'
    };
    % Parietal/Temporal Nodes (P/T) - Represents Parietal/Temporal Lobe Activity
    posterior_nodes = {
        'TP7-CP5', 'CP5-P7', 'P7-PO7', 'Pz-POz', 'P3-P5', 'P4-P6', 'CP1-P3', 'CP2-P4', 'C6-TP8', 'TP8-CP6', 'CP6-P8', 'P8-PO8', 'PO8-O2'
    };
end

%%
%% ========= CONNECTIVITY ANALYSIS BLOCK (TARGETED F-P NETWORK) ============
fprintf('\n=== OPTIMIZED CONNECTIVITY ANALYSIS (TARGETED F-P NETWORK PLV) ===\n');

% Get the full list of bipolar channel labels
bipolar_labels = {CONN_EEG.chanlocs.labels};
[frontal_node_labels, posterior_node_labels] = get_target_network_groups();

% --- CRITICAL STEP 1: Map Node Labels to Indices in CONN_EEG ---
frontal_indices = [];
for i = 1:length(frontal_node_labels)
    idx = find(strcmpi(bipolar_labels, frontal_node_labels{i}));
    if ~isempty(idx)
        frontal_indices(end+1) = idx;
    else
        fprintf('Warning: Frontal node %s not found. Skipping.\n', frontal_node_labels{i});
    end
end

posterior_indices = [];
for i = 1:length(posterior_node_labels)
    idx = find(strcmpi(bipolar_labels, posterior_node_labels{i}));
    if ~isempty(idx)
        posterior_indices(end+1) = idx;
    else
        fprintf('Warning: Posterior node %s not found. Skipping.\n', posterior_node_labels{i});
    end
end

if isempty(frontal_indices) || isempty(posterior_indices)
    error('Insufficient frontal or posterior nodes found for F-P connectivity analysis.');
end

% The selected channels are all channels that belong to either F or P/T nodes.
all_target_indices = unique([frontal_indices, posterior_indices]);

fprintf('Targeted analysis using %d frontal nodes and %d posterior nodes.\n', ...
    length(frontal_indices), length(posterior_indices));


%% 1. DEFINE CONDITIONS
conditions = {
    'startOfRememberedClipFirstWatch', 'Remembered';
    'startOfRecognisedClipFirstWatch', 'Recognised';
    'startOfNotRecognisedClip', 'Forgotten'
};
conditions_to_test = {'Remembered','Recognised','Forgotten'};

%% 2. Extract data for each condition by EPOCH INDEX (and Subset Channels)
fprintf('Preparing condition datasets and extracting data by epoch index...\n');
condition_data = struct();
epoch_types_all = cell(1, length(CONN_EEG.epoch));
for i = 1:length(CONN_EEG.epoch)
    if iscell(CONN_EEG.epoch(i).eventtype), epoch_types_all{i} = CONN_EEG.epoch(i).eventtype{1};
    else, epoch_types_all{i} = CONN_EEG.epoch(i).eventtype; end
end

for i = 1:size(conditions,1)
    cond_name  = conditions{i,2};
    event_name = conditions{i,1};
    trial_indices = find(strcmpi(epoch_types_all, event_name));
    
    if ~isempty(trial_indices)
        temp_EEG = pop_select(CONN_EEG, 'trial', trial_indices);
        
        % CRITICAL STEP 2: Select ONLY the target channels for the PLV calculation
        temp_EEG = pop_select(temp_EEG, 'channel', all_target_indices);
        
        % CRITICAL STEP 3: Find the indices of the F and P nodes *within this new subset EEG*
        % Since pop_select preserves the order of the indices given in all_target_indices,
        % we need to know the mapping of the F and P nodes within the target subset.
        
        % A more robust way: use the original bipolar labels and match them against the subset labels.
        subset_labels = {temp_EEG.chanlocs.labels};
        
        subset_frontal_indices = [];
        for j = 1:length(frontal_node_labels)
            idx = find(strcmpi(subset_labels, frontal_node_labels{j}));
            if ~isempty(idx)
                subset_frontal_indices(end+1) = idx;
            end
        end

        subset_posterior_indices = [];
        for j = 1:length(posterior_node_labels)
            idx = find(strcmpi(subset_labels, posterior_node_labels{j}));
            if ~isempty(idx)
                subset_posterior_indices(end+1) = idx;
            end
        end
        
        condition_data.(cond_name).data = temp_EEG.data;
        condition_data.(cond_name).trials = temp_EEG.trials;
        condition_data.(cond_name).F_indices = subset_frontal_indices;
        condition_data.(cond_name).P_indices = subset_posterior_indices;
        
        fprintf('%s: %d trials extracted (%d total channels in F+P subset).\n', cond_name, temp_EEG.trials, temp_EEG.nbchan);
    else
        condition_data.(cond_name) = [];
        fprintf('%s: NO trials found.\n', cond_name);
    end
end
% --- End of Extraction ---

%% 3. Calculate minimum trials and Subsample (Equalize)
extracted_trials = zeros(1, length(conditions_to_test));
for c_idx = 1:length(conditions_to_test)
    cond = conditions_to_test{c_idx};
    if isfield(condition_data, cond) && ~isempty(condition_data.(cond))
        extracted_trials(c_idx) = condition_data.(cond).trials;
    else
        extracted_trials(c_idx) = 0;
    end
end
minTrials = min(extracted_trials(extracted_trials > 0));

if isempty(minTrials) || minTrials == 0
    warning('Connectivity: Minimum trial count is 0. Skipping PLV calculation.');
    PART_CONN_SKIP = true;
else
    PART_CONN_SKIP = false;    
    
    fprintf('Equalizing trials to minimum count: %d\n', minTrials);
    for i = 1:length(conditions_to_test)
        cond = conditions_to_test{i};
        if ~isempty(condition_data.(cond)) && condition_data.(cond).trials > minTrials
            idx = randperm(condition_data.(cond).trials, minTrials);
            condition_data.(cond).data = condition_data.(cond).data(:,:,idx);
            condition_data.(cond).trials = minTrials;
            fprintf('%s: reduced to %d trials.\n', cond, minTrials);
        end
    end

    %% 4. PRE-CALCULATE PLV PER TRIAL (OPTIMIZED)
    fprintf('\n=== STEP 4: PRE-CALCULATING PLV FOR EACH TRIAL (4-8 Hz) IN F-P NETWORK ===\n');
    theta_trial_plvs = struct();
    theta_results = struct();
    for i = 1:length(conditions_to_test)
        cond = conditions_to_test{i};
        if ~isempty(condition_data.(cond))
            EEGtmp = struct('data', condition_data.(cond).data, 'srate', CONN_EEG.srate);
            F_idx = condition_data.(cond).F_indices;
            P_idx = condition_data.(cond).P_indices;
            
            % PLV calculation now only runs between F_idx and P_idx
            plvs = plv_per_trial(EEGtmp, [4 8], F_idx, P_idx); 
            
            theta_trial_plvs.(cond) = plvs;
            theta_results.(cond).avg_plv = mean(plvs);
            theta_results.(cond).ntrials = condition_data.(cond).trials;
            
            num_fp_pairs = length(F_idx) * length(P_idx);
            fprintf('%s: AVG F-P PLV = %.3f (n=%d, %d F-P pairs)\n', cond, theta_results.(cond).avg_plv, theta_results.(cond).ntrials, num_fp_pairs);
        else
            theta_results.(cond).avg_plv = NaN;
            theta_results.(cond).ntrials = 0;
        end
    end

    %% 5. FAST PERMUTATION TEST BETWEEN CONDITIONS
    fprintf('\n=== STEP 5: FAST PERMUTATION TEST (N=1000) ===\n');
    perm_pairs = {'Remembered','Forgotten'; 'Remembered','Recognised'; 'Recognised','Forgotten'};
    nPerm = 1000; 
    perm_results = struct();
    for i = 1:size(perm_pairs,1)
        c1 = perm_pairs{i,1}; c2 = perm_pairs{i,2};
        
        if ~isfield(theta_trial_plvs, c1) || ~isfield(theta_trial_plvs, c2) || isempty(theta_trial_plvs.(c1)) || isempty(theta_trial_plvs.(c2))
            fprintf('%s vs %s: missing trial PLVs, skipping.\n', c1, c2);
            continue
        end
        
        plv1 = theta_trial_plvs.(c1);
        plv2 = theta_trial_plvs.(c2);
        n1 = length(plv1); n2 = length(plv2);
        
        if n1 ~= n2 || n1 == 0 
            fprintf('%s vs %s: zero or unequal trials after subsampling (%d vs %d), skipping.\n', c1, c2, n1, n2);
            continue
        end
        
        obs_diff = theta_results.(c1).avg_plv - theta_results.(c2).avg_plv;
        combined_plvs = [plv1, plv2];
        perm_diff = zeros(1, nPerm);
        
        for p = 1:nPerm
            idx = randperm(n1+n2);
            perm_group1_plv = combined_plvs(idx(1:n1));
            perm_group2_plv = combined_plvs(idx(n1+1:end));
            perm_diff(p) = mean(perm_group1_plv) - mean(perm_group2_plv);
        end
        
        p_val = mean(abs(perm_diff) >= abs(obs_diff));
        fprintf('%s vs %s: observed diff = %.3f, p = %.4f\n', c1, c2, obs_diff, p_val);
        perm_results.(sprintf('%s_vs_%s', c1, c2)).obs_diff = obs_diff;
        perm_results.(sprintf('%s_vs_%s', c1, c2)).p_val = p_val;
        perm_results.(sprintf('%s_vs_%s', c1, c2)).perm_diff = perm_diff;
    end
end
% --- End of Connectivity Analysis Block ---

%% ========= SUPPORTING FUNCTIONS (Optimization Helpers) ================
function plv_trial = plv_per_trial(EEG, band, frontal_indices, posterior_indices)
    % Calculates the average PLV ONLY between Frontal and Posterior nodes.
    data = EEG.data; 
    [nch, ~, ntrials] = size(data);
    srate = EEG.srate;
    
    filtered = bandpass_filter(data, srate, band);
    plv_trial = zeros(1, ntrials);
    
    num_f_nodes = length(frontal_indices);
    num_p_nodes = length(posterior_indices);
    total_fp_pairs = num_f_nodes * num_p_nodes;
    
    if total_fp_pairs == 0
        plv_trial = zeros(1, ntrials);
        return;
    end

    for tr = 1:ntrials
        phase_tr = angle(hilbert(squeeze(filtered(:,:,tr))'))';
        complex_phase_tr = exp(1i * phase_tr);
        
        plv_sum = 0;
        
        % Compute PLV ONLY between F nodes and P nodes
        for ch1_idx = 1:num_f_nodes
            ch1 = frontal_indices(ch1_idx); % Index into the data matrix
            for ch2_idx = 1:num_p_nodes
                ch2 = posterior_indices(ch2_idx); % Index into the data matrix
                
                plv_sum = plv_sum + abs(mean(complex_phase_tr(ch1, :) .* conj(complex_phase_tr(ch2, :))));
            end
        end
        plv_trial(tr) = plv_sum / total_fp_pairs;
    end
end

function out = bandpass_filter(data, srate, band)
    % 4th order Butterworth bandpass filter applied to 3D EEG data
    [b,a] = butter(4, band/(srate/2), 'bandpass');
    
    out = zeros(size(data));
    for ch = 1:size(data, 1)
        for tr = 1:size(data, 3)
            out(ch,:,tr) = filtfilt(b, a, double(squeeze(data(ch,:,tr))));
        end
    end
end
fprintf('\n=== FULL PIPELINE EXECUTION COMPLETE ===\n');