% Step 1: Load sample EEG data
load sampleEEGdata;   % loads EEG structure

% Select one electrode (different from example, e.g., 'Pz')
chanName = 'Pz';
chanIdx = find(strcmpi({EEG.chanlocs.labels}, chanName));

% Extract data and parameters
data = EEG.data(chanIdx, :, :);   % channel x time x trials
fs = EEG.srate;                   % sampling rate
times = EEG.times;                % time vector

% Average across trials to simplify visualization
signal = mean(data, 3);

%% Complex morelet wavelet convolution
% Step 2: Power at 10 Hz using Complex Morlet Wavelet
freq = 10;              % Frequency of interest (Hz)
fs = EEG.srate;

t = -2:1/fs:2;          % 4-second window centered at 0
s = 4 / (2*pi*freq);    % Gaussian width (4 cycles)
wavelet = exp(2i*pi*freq*t) .* exp(-t.^2/(2*s^2));

convResult = conv(signal, wavelet, 'same');
power_wavelet = abs(convResult).^2;

figure;
plot(times, power_wavelet, 'LineWidth', 1.5)
xlabel('Time (ms)')
ylabel('Power (a.u.)')
title(['Power at ', num2str(freq), ' Hz — Complex Wavelet'])

%% Filter-Hilbert Method
% Step 3: Power at 10 Hz using Filter–Hilbert method
freq = 10;
band = [8 12];   % narrow band around 10 Hz
fs = EEG.srate;

% Band-pass filter (2nd order Butterworth)
[b,a] = butter(2, band/(fs/2), 'bandpass');
filtSignal = filtfilt(b, a, double(signal));

% Hilbert transform to get analytic signal
analytic = hilbert(filtSignal);
power_hilbert = abs(analytic).^2;

figure;
plot(times, power_hilbert, 'r', 'LineWidth', 1.5)
xlabel('Time (ms)')
ylabel('Power (a.u.)')
title(['Power at ', num2str(freq), ' Hz — Filter–Hilbert'])

%% STFT Method
% Step 4: Power at 10 Hz using Short-Time FFT
win_len = round(0.5 * fs);     % 500 ms window
step    = round(0.05 * fs);    % 50 ms step
nfft    = 512;
freq    = 10;

times_stft = [];
power_stft = [];

for start = 1:step:(length(signal)-win_len)
    segment = signal(start:start+win_len-1);
    Y = fft(segment, nfft);
    f = (0:nfft-1)*(fs/nfft);
    power_spectrum = abs(Y).^2 / nfft;
    [~, idx] = min(abs(f - freq));
    power_stft(end+1) = power_spectrum(idx);
    times_stft(end+1) = times(start + win_len/2);
end

figure;
plot(times_stft, power_stft, 'b', 'LineWidth', 1.5)
xlabel('Time (ms)')
ylabel('Power (a.u.)')
title(['Power at ', num2str(freq), ' Hz — Short-Time FFT'])




%% Multi taper Method
% Step 5: Power at 10 Hz using Multitaper method

win_len = round(0.5 * fs);     % 500 ms window
step    = round(0.05 * fs);    % 50 ms step
nw      = 3;                   % Time-halfbandwidth product
nfft    = 512;
freq    = 10;

times_mt = [];
power_mt = [];

for start = 1:step:(length(signal)-win_len)
    segment = signal(start:start+win_len-1);
    [Pxx, f] = pmtm(segment, nw, nfft, fs);
    [~, idx] = min(abs(f - freq));
    power_mt(end+1) = Pxx(idx);
    times_mt(end+1) = times(start + win_len/2);
end

figure;
plot(times_mt, power_mt, 'm', 'LineWidth', 1.5)
xlabel('Time (ms)')
ylabel('Power (a.u.)')
title(['Power at ', num2str(freq), ' Hz — Multitaper'])

%% Baseline Comparision and Comaprision 
% Step 6: Baseline correction and comparison of all methods

% Define baseline range (pre-stimulus)
base_idx = times >= -500 & times <= 0;

% Compute mean baseline power for each method
base_wavelet = mean(power_wavelet(base_idx));
base_hilbert = mean(power_hilbert(base_idx));
base_stft    = mean(power_stft(times_stft >= -500 & times_stft <= 0));
base_mt      = mean(power_mt(times_mt >= -500 & times_mt <= 0));

% Normalize (divide baseline)
wavelet_norm = power_wavelet / base_wavelet;
hilbert_norm = power_hilbert / base_hilbert;
stft_norm    = power_stft / base_stft;
mt_norm      = power_mt / base_mt;

% Plot all normalized power curves
figure; hold on
plot(times, wavelet_norm, 'k', 'LineWidth', 1.5)
plot(times, hilbert_norm, 'r', 'LineWidth', 1.5)
plot(times_stft, stft_norm, 'b', 'LineWidth', 1.5)
plot(times_mt, mt_norm, 'm', 'LineWidth', 1.5)
xlabel('Time (ms)')
ylabel('Normalized Power')
title('10 Hz Power — All Methods (Baseline-corrected)')
legend('Wavelet', 'Hilbert', 'STFT', 'Multitaper')
