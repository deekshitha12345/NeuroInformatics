% eeglab's EEG structure and indexing
% mikeXcohen@gmail.com

%% load in EEG data

load sampleEEGdata.mat

% FYI, this would also work:
% [file2load,path4file]=uigetfile('*.mat','Please select EEG data file');
% load([ path4file file2load ])

% take a minute to inspect the EEG structure
EEG
whos EEG                 % show EEG size and type
fieldnames(EEG)          % list the fields in the struct
size(EEG.data)           % should return [nbchan, pnts, trials]
EEG.nbchan               % number of channels
EEG.pnts                 % number of time points (per trial)
EEG.trials               % number of trials / epochs
EEG.srate                % sampling rate (Hz)
EEG.times(1:10)          % first 10 time samples (units usually ms)
{EEG.chanlocs.labels}    % cell array of channel names
openvar("EEG")           % To open struct
 

%% Finding time indices based on ms

time2plot = 300; % desired time in ms

% find closest index in EEG.times
[~,timeidx] = min(abs(EEG.times - time2plot));

% extract the data at that time across trials
% (here trial average across 3rd dim)
data2plot = squeeze(mean(EEG.data(:,timeidx,:) ,3));

% Plot with topoplotIndie (instead of topoplot)
figure(1), clf
topoplotIndie(data2plot, EEG.chanlocs);

% Add title showing both the requested and actual time
title([ 'Topoplot at requested ' num2str(time2plot) ' ms (closest = ' num2str(EEG.times(timeidx)) ' ms)' ])


%% same concept for frequencies

frex = linspace(2,100,42);

freqIwant = 23; % in hz

% use min(abs trick to find closest frequency to 23 Hz
[~,frexidx] = min(abs(frex-freqIwant));

% the function dsearchn also works
frexidx = dsearchn(frex',freqIwant);

%% indexing channels based on names

% the electrode label that we want to analyze
electrodeName = 'p1'; % case doesn't matter

% find the channel number that corresponds to this label
electrodeidx = strcmpi(electrodeName,{EEG.chanlocs.labels});

% confirm that this electrode is correct
EEG.chanlocs(electrodeidx)

% plot the ERP from this electrode
figure(1), clf
plot(EEG.times,mean( EEG.data(electrodeidx,:,:),3 ))


%% now multiple electrodes


electrodeNames = {'p1','fc6','t18'};

% initialize
electrodeidx = zeros(1,length(electrodeNames));

% loop through electrodes and find the index of each one

for chani=1:length(electrodeNames)
    try
        % strcmpi gives a logical array → convert it to index using find()
        outputarray = strcmpi(electrodeNames{chani},{EEG.chanlocs.labels});
        electrodeidx(chani) = find(outputarray);
    catch ME
        % if electrode is not found, store NaN and print warning
        disp(['Electrode "' electrodeNames{chani} '" not found.']);
        electrodeidx(chani) = NaN;
    end
end

% remove NaNs before indexing into EEG.data
electrodeidx(isnan(electrodeidx)) = [];

% if valid electrodes found, plot ERP
if ~isempty(electrodeidx)
    figure, plot(EEG.times, mean(EEG.data(electrodeidx,:,:),3))
    legend({EEG.chanlocs(electrodeidx).labels})
    xlabel('Time (ms)'), ylabel('Amplitude (\muV)')
    title('ERP from selected electrodes')
else
    disp('No valid electrodes found. Nothing to plot.');
end


%%

