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

%% Step 1: Downsampling
fprintf('\n=== STEP 1: DOWNSAMPLING ===\n');
% Explanation: Reduces sampling rate from 2048 Hz to 512 Hz
% - Speeds up processing (4x fewer data points)
% - Still captures all relevant brain frequencies (< 200 Hz)
% - Reduces computational load for ICA
original_srate = EEG.srate;
EEG = pop_resample(EEG, 512);
EEG.setname = 'Downsampled_512Hz';
EEG = eeg_checkset(EEG);
fprintf('Downsampling complete: %d Hz → %d Hz\n', original_srate, EEG.srate);

%% Step 2: High Pass Filtering
fprintf('\n=== STEP 2: HIGH PASS FILTERING ===\n');
% Explanation: Removes slow drifts and DC offsets
% - 0.1 Hz cutoff removes very slow changes
% - Essential for ICA - it needs stationary data
% - Removes baseline shifts that could confuse ICA
EEG = pop_eegfiltnew(EEG, 0.1, [], [], false, [], 0);
EEG.setname = 'High_Pass_Filtered';
EEG = eeg_checkset(EEG);
fprintf('High-pass filtering complete: 0.1 Hz cutoff\n');

%% Step 3: Bad Channel Detection and Removal
fprintf('\n=== STEP 3: BAD CHANNEL DETECTION ===\n');
% Explanation: Removes channels with poor signal quality
% - Flat channels (no activity for >5 seconds)
% - Channels poorly correlated with neighbors (<0.8)
% - Channels with excessive 50/60 Hz noise
% - Bad channels would create bad ICA components
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
if ~isempty(removed_channels)
    fprintf('  Removed channels: %s\n', strjoin(removed_channels, ', '));
else
    fprintf('  No channels removed\n');
end

%% Step 4: Average Referencing  
fprintf('\n=== STEP 4: AVERAGE REFERENCING ===\n');
% Explanation: Changes reference from unknown to average of all channels
% - Removes common noise across all electrodes
% - Improves ICA performance by removing shared signals
% - Standard preprocessing step for ICA
EEG = pop_reref(EEG, []);
EEG.setname = 'Average_Referenced';
EEG = eeg_checkset(EEG);
fprintf('Average referencing complete\n');

%% Step 5: Low Pass Filtering
fprintf('\n=== STEP 5: LOW PASS FILTERING ===\n');
% Explanation: Removes high-frequency noise
% - 40 Hz cutoff removes muscle artifacts and electrical noise
% - Keeps all relevant brain signals (most are < 40 Hz)
% - Clean data improves ICA decomposition quality
EEG = pop_eegfiltnew(EEG, [], 40, [], false, [], 0);
EEG.setname = 'Band_Pass_Filtered';
EEG = eeg_checkset(EEG);
fprintf('Low-pass filtering complete: 40 Hz cutoff\n');

%% Step 6: ICA PREPARATION AND COMPUTATION
fprintf('\n=== STEP 6: ICA DECOMPOSITION ===\n');

% Save pre-ICA data for comparison
EEG_before_ICA = EEG;

% IMPORTANT: ICA needs continuous data, so we do it BEFORE epoching
fprintf('Running ICA decomposition...\n');
fprintf('This may take several minutes for %d channels...\n', EEG.nbchan);

% Run ICA using pop_runica (EEGLAB's ICA function)
% 'runica' is the default algorithm (Infomax ICA)
% 'extended', 1 handles sub-gaussian and super-gaussian sources
EEG = pop_runica(EEG, 'icatype', 'runica', 'extended', 1);
EEG.setname = 'ICA_Complete';
EEG = eeg_checkset(EEG);

fprintf('ICA decomposition complete!\n');
fprintf('Number of ICA components: %d\n', size(EEG.icawinv, 2));

%% Step 7: VISUALIZE ALL ICA COMPONENTS
fprintf('\n=== STEP 7: ICA COMPONENT ANALYSIS ===\n');

% Get total number of components
n_comps = size(EEG.icawinv, 2);
fprintf('Total ICA components: %d\n', n_comps);

% Show ALL components (not just 35)
try
    fprintf('Displaying ALL %d component topographies...\n', n_comps);
    pop_selectcomps(EEG, 1:n_comps);  % Show ALL components
    fprintf('ICA component topographies displayed successfully.\n');
catch ME
    fprintf('Error displaying component topographies: %s\n', ME.message);
    fprintf('Trying alternative visualization...\n');
    
    % Alternative: Manual topography plotting
    try
        fprintf('Creating manual component visualization...\n');
        
        % Calculate grid size for subplots
        n_rows = ceil(sqrt(n_comps));
        n_cols = ceil(n_comps / n_rows);
        
        figure('Name', 'ICA Component Topographies', 'Position', [50 50 1400 1000]);
        for comp_i = 1:n_comps
            subplot(n_rows, n_cols, comp_i);
            
            % Use topoplotIndie if available, otherwise skip topography
            try
                if exist('topoplotIndie', 'file')
                    topoplotIndie(EEG.icawinv(:, comp_i), EEG.chanlocs, ...
                        'electrodes', 'off', 'shading', 'interp');
                else
                    % Simple bar plot of component weights
                    bar(EEG.icawinv(:, comp_i));
                    xlabel('Channel');
                end
                title(sprintf('IC %d', comp_i));
            catch
                text(0.5, 0.5, sprintf('IC %d\n(plot error)', comp_i), ...
                    'HorizontalAlignment', 'center');
                axis off;
            end
        end
        sgtitle('All ICA Components');
        
    catch ME2
        fprintf('Could not create manual visualization: %s\n', ME2.message);
        fprintf('ICA was computed successfully, but visualization failed.\n');
    end
end

% Show component time courses
try
    fprintf('Displaying component time courses...\n');
    pop_eegplot(EEG, 0, 1, 1); % 0 = plot components
    fprintf('Component time courses displayed.\n');
catch ME
    fprintf('Could not display component time courses: %s\n', ME.message);
end

fprintf('\nLook for these artifact patterns:\n');
fprintf('  - Eye blinks: Strong frontal (Fp1, Fp2) activation\n');
fprintf('  - Eye movements: Asymmetric frontal activation\n');
fprintf('  - Muscle: Localized temporal activation\n');
fprintf('  - Heart: Lower channels, rhythmic pattern\n');
fprintf('  - Line noise: Distributed, 50/60 Hz frequency\n');

%% Step 8: MANUAL COMPONENT SELECTION
fprintf('\n=== STEP 8: MANUAL COMPONENT SELECTION ===\n');

% Display component properties for better identification
try
    fprintf('Computing component properties...\n');
    EEG = pop_companalyze(EEG);
    fprintf('Component analysis complete.\n');
catch
    fprintf('Component analysis not available, proceeding with basic selection.\n');
end

% Enhanced component information
fprintf('\nComponent Information:\n');
fprintf('Total components: %d\n', n_comps);
fprintf('Components to inspect: 1 to %d\n', n_comps);

% Show component selection interface again
try
    pop_selectcomps(EEG, 1:n_comps);  % Show ALL components for selection
catch
    fprintf('GUI component selection failed.\n');
end

% Manual selection with validation
fprintf('\n=== COMPONENT REMOVAL SELECTION ===\n');
fprintf('Based on the component topographies and time courses:\n');
fprintf('1. Identify components that look like artifacts\n');
fprintf('2. Eye blinks: usually IC1-IC3, strong frontal\n');
fprintf('3. Muscle: localized, high frequency\n');
fprintf('4. Heart: rhythmic, lower electrodes\n');
fprintf('5. Line noise: 60Hz, distributed\n\n');

% Get user input with validation
valid_input = false;
while ~valid_input
    fprintf('Enter component numbers to remove:\n');
    fprintf('Examples: [1 2 5] or [1:3 7] or [] for none\n');
    fprintf('Your choice: ');
    
    try
        user_input = input('');
        
        if isempty(user_input)
            manual_components = [];
            valid_input = true;
            fprintf('No components selected for removal.\n');
        elseif isnumeric(user_input) && all(user_input >= 1) && all(user_input <= n_comps)
            manual_components = unique(user_input);
            valid_input = true;
            fprintf('Selected components: %s\n', mat2str(manual_components));
        else
            fprintf('Invalid input. Please enter numbers between 1 and %d.\n', n_comps);
        end
    catch
        fprintf('Invalid input format. Please try again.\n');
    end
end

% Remove selected components
if ~isempty(manual_components)
    fprintf('\nRemoving %d components: %s\n', length(manual_components), mat2str(manual_components));
    
    % Show what will be removed
    for i = 1:length(manual_components)
        comp_num = manual_components(i);
        fprintf('  Removing IC%d\n', comp_num);
    end
    
    % Confirm removal
    fprintf('\nProceed with removal? (y/n): ');
    confirm = input('', 's');
    
    if strcmpi(confirm, 'y')
        EEG = pop_subcomp(EEG, manual_components, 0);
        fprintf('Successfully removed %d artifact components.\n', length(manual_components));
        all_artifact_components = manual_components;
    else
        fprintf('Component removal cancelled.\n');
        all_artifact_components = [];
    end
else
    fprintf('No components removed.\n');
    all_artifact_components = [];
end

%% Step 9: MANUAL COMPONENT SELECTION (if needed)
if ~exist('all_artifact_components', 'var') || isempty(all_artifact_components)
    fprintf('\n=== STEP 9: MANUAL COMPONENT SELECTION ===\n');
    
    % Display component selection interface
    fprintf('Please manually select artifact components:\n');
    fprintf('1. Look at the component topographies\n');
    fprintf('2. Identify artifact patterns (eye, muscle, heart, noise)\n');
    fprintf('3. Enter component numbers to remove\n');
    
    % Show component selection GUI again if needed
    pop_selectcomps(EEG, 1:min(35, size(EEG.icawinv, 2)));
    
    % Manual input for component removal
    fprintf('\nEnter component numbers to remove (e.g., [1 3 5] or [] for none): ');
    manual_components = input('');
    
    if ~isempty(manual_components)
        EEG = pop_subcomp(EEG, manual_components, 0);
        all_artifact_components = manual_components;
        fprintf('Removed %d components manually: %s\n', ...
            length(manual_components), mat2str(manual_components));
    else
        fprintf('No components removed manually.\n');
        all_artifact_components = [];
    end
end

%% Step 10: COMPARE BEFORE AND AFTER ICA
fprintf('\n=== STEP 10: ICA RESULTS COMPARISON ===\n');

% Create comparison figure
figure('Name', 'ICA Results Comparison', 'Position', [100 100 1400 800]);

% Plot a few channels before and after ICA
channels_to_plot = [1, 8, 16, 32]; % Sample channels
time_window = [1, min(30*EEG.srate, EEG.pnts)]; % First 30 seconds
time_axis = (time_window(1):time_window(2)) / EEG.srate;

for i = 1:length(channels_to_plot)
    ch = channels_to_plot(i);
    
    % Before ICA
    subplot(length(channels_to_plot), 2, 2*i-1);
    plot(time_axis, EEG_before_ICA.data(ch, time_window(1):time_window(2)));
    title(sprintf('Channel %s - Before ICA', EEG.chanlocs(ch).labels));
    ylabel('μV');
    if i == 1
        title(sprintf('Channel %s - Before ICA (with artifacts)', EEG.chanlocs(ch).labels));
    end
    
    % After ICA
    subplot(length(channels_to_plot), 2, 2*i);
    plot(time_axis, EEG.data(ch, time_window(1):time_window(2)));
    title(sprintf('Channel %s - After ICA', EEG.chanlocs(ch).labels));
    ylabel('μV');
    if i == 1
        title(sprintf('Channel %s - After ICA (artifacts removed)', EEG.chanlocs(ch).labels));
    end
    
    if i == length(channels_to_plot)
        subplot(length(channels_to_plot), 2, 2*i-1);
        xlabel('Time (s)');
        subplot(length(channels_to_plot), 2, 2*i);
        xlabel('Time (s)');
    end
end

sgtitle('ICA Artifact Removal Results');

%% Step 11: Epoching
fprintf('\n=== STEP 11: EPOCHING ===\n');
% Now we can epoch the clean data

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

%% Step 12: Baseline Correction
fprintf('\n=== STEP 12: BASELINE CORRECTION ===\n');
% Remove baseline using pre-stimulus period
EEG = pop_rmbase(EEG, [-200 0]);
EEG.setname = 'Baseline_Corrected';
EEG = eeg_checkset(EEG);
fprintf('Baseline correction complete using [-200, 0] ms\n');

%% Step 13: Final Artifact Rejection
fprintf('\n=== STEP 13: FINAL ARTIFACT REJECTION ===\n');
% Remove any remaining bad epochs
epochs_before = EEG.trials;

EEG = pop_autorej(EEG, 'nogui', 'on', 'threshold', 100, 'startprob', 5, ...
                  'maxrej', 5, 'eegplot', 'off');

epochs_after = EEG.trials;
epochs_rejected = epochs_before - epochs_after;

fprintf('Final artifact rejection complete:\n');
fprintf('  Epochs before: %d\n', epochs_before);
fprintf('  Epochs after: %d\n', epochs_after);
fprintf('  Epochs rejected: %d (%.1f%%)\n', epochs_rejected, ...
    (epochs_rejected/epochs_before)*100);

%% Final Summary with ICA Information
fprintf('\n=== PREPROCESSING SUMMARY WITH ICA ===\n');
fprintf('Original sampling rate: 2048 Hz\n');
fprintf('Final sampling rate: %d Hz\n', EEG.srate);
fprintf('Original channels: 64\n');
fprintf('Channels after cleaning: %d\n', EEG.nbchan);
if ~isempty(removed_channels)
    fprintf('Removed channels: %s\n', strjoin(removed_channels, ', '));
else
    fprintf('No bad channels removed\n');
end
fprintf('ICA components computed: %d\n', size(EEG.icawinv, 2));
if exist('all_artifact_components', 'var') && ~isempty(all_artifact_components)
    fprintf('ICA components removed: %d (%.1f%%)\n', length(all_artifact_components), ...
        length(all_artifact_components)/size(EEG.icawinv, 2)*100);
    fprintf('Removed component numbers: %s\n', mat2str(all_artifact_components));
else
    fprintf('No ICA components removed\n');
end
fprintf('Total epochs: %d\n', EEG.trials);

%% Save processed data with ICA information
save_filename = sprintf('%s_task-%s_preprocessed_withICA.set', subject, task);
pop_saveset(EEG, 'filename', save_filename, 'filepath', pwd);
fprintf('\nPreprocessed data with ICA saved as: %s\n', save_filename);

fprintf('\n=== PREPROCESSING WITH ICA COMPLETE ===\n');
fprintf('Your data is now clean and ready for analysis!\n');



