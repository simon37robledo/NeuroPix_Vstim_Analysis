% Use Welch’s method to estimate power spectral density
[pxx, f] = pwelch(Tempdata, [], [], [], 99);

% Plot
figure;
plot(f,10*log10(pxx));
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density of DLC coordinate');

%%
fs = 100;             % Sampling frequency (Hz)
signal = Tempdata;      % Example coordinate

% Parameters for Welch's method
windowLength = 10*fs; % 10-second window = 1000 samples
noverlap = floor(windowLength/2);
nfft = 8192;          % Higher nfft = better frequency resolution

[pxx, f] = pwelch(signal, windowLength, noverlap, nfft, fs);

% Plot
semilogx(f, 10*log10(pxx));  % Use log scale to see low freqs better
xlim([0.01 5]);              % Focus on <5 Hz
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('PSD showing low-frequency fluctuations');

