%load only the file you want to inspect

%load phantomRFdata1;
%p
%load phantomRFdata2;
%p
load cardiacRfdata;
p


% Task 1

% 'fs' (sampling frequency) and 'f0' (transmit center frequency)
fs = p.frs_Hz; 
f0 = p.f0_Hz;
rfData = rf;

% Design a low-pass filter
lowPassFilter = fir1(60, 0.2, "low");



% Call the complexDemodulation function
IQdata = complexDemodulation(rfData, f0, fs, lowPassFilter);

harmonicIQ = complexDemodulationHarmonic(rfData, f0, fs); %only relevant for cardiac

%Calculate power spectra
Nfft = 2048;
PowerSpectrum_RF=20*log10(abs(fftshift(fft(rfData, Nfft)))); %Power spectrum unmodified
AveragePowerSpectrum_RF = mean(PowerSpectrum_RF, 2); % Average over all columns

PowerSpectrum_IQ=20*log10(abs(fftshift(fft(IQdata, Nfft)))); %Power spectrum demodulated
AveragePowerSpectrum_IQ = mean(PowerSpectrum_IQ, 2); % Average over all columns

PowerSpectrum_harmonicIQ=20*log10(abs(fftshift(fft(harmonicIQ, Nfft)))); %Power spectrum demodulated cardiac
AveragePowerSpectrum_harmonicIQ = mean(PowerSpectrum_harmonicIQ, 2); % Average over all columns cardiac

% Corrected frequency vector for plotting
f = linspace(-fs/2, fs/2, Nfft);

% Adjusted plotting to use the corrected frequency vector
% And ensure plotting across the correct frequency range
figure;
plot(f, AveragePowerSpectrum_RF, 'b', 'LineWidth', 2); % Original RF data
hold on;
plot(f, AveragePowerSpectrum_IQ, 'r', 'LineWidth', 2); % Demodulated IQ data
hold on;
plot(f, AveragePowerSpectrum_harmonicIQ, 'g', 'LineWidth', 2); % Demodulated IQ data only relevant for cardiac
xlabel('Frequency (Hz)');
ylabel('Power Spectrum (dB)');
legend('Original RF', 'Demodulated IQ', 'Demodulated IQ Harmonic');
title('Power Spectrum: Original RF vs. Demodulated IQ vs. Demodulated IQ Harmonic');
grid on;



% Task 2
c = 1540;

% Calculate the middle image line index
middleIndex = ceil(size(rfData, 2) / 2);

% Extract the RF signal from the middle image line
middleRFSignal = rfData(:, middleIndex);

% Calculate the depth axis
% Time for each sample (in seconds)
time = (0:size(rfData, 1) - 1)'/fs;
% Depth calculation
depth = c * time / 2 + p.startdepth_m;

% Plot the RF signal as a function of depth
figure;
plot(depth, middleRFSignal);
xlabel('Depth (m)');
ylabel('RF Amplitude');
title('RF Signal from Middle Image Line as a Function of Depth');
grid on;

% Envelope
envelopeIQ = 2 * abs(IQdata); % Calculate envelope and multiply by 2

% Extract the envelope of the middle image line
middleEnvelope = envelopeIQ(:, middleIndex);

% Plotting
figure;
plot(depth, middleRFSignal); % Original RF signal from 2a
hold on;
plot(depth, middleEnvelope, 'r', 'LineWidth', 0.1); % Envelope of the middle image line
xlabel('Depth (m)');
ylabel('Signal Amplitude');
legend('Original RF Signal', 'Envelope');
title('RF Signal and Its Envelope from the Middle Image Line');
grid on;



% Task 3
num_beams = size(IQdata, 2);

% Calculate theta vector
theta = p.startangle_rad + (0:(num_beams-1)) * p.angleincrement_rad;

% Calculate time for each sample
time_vector = (0:(size(IQdata, 1) - 1))' / fs;

% Calculate depth vector 'r', including start depth
r = c * time_vector / 2 + p.startdepth_m;

% Perform scan conversion

[cartesianImage, x_axis, y_axis] = scanconvert(IQdata, r, theta);

logCompressedImage = 20 * log10(abs(cartesianImage));

% Define dynamic range and gain
dyn = 60; % Dynamic range in dB
gain = -100; 

% Display the log-compressed image
figure;
imagesc(logCompressedImage);
caxis([-dyn 0] - gain); % Adjust contrast and brightness
colormap(gray); % Set colormap to gray for ultrasound images
colorbar;

xlabel('Lateral Distance');
ylabel('Depth');
title('Log-Compressed Ultrasound Image');


%Task 5
num_beams = size(harmonicIQ, 2);

% Calculate theta vector
theta = p.startangle_rad + (0:(num_beams-1)) * p.angleincrement_rad;

% Calculate time for each sample
time_vector = (0:(size(harmonicIQ, 1) - 1))' / fs;

% Calculate depth vector 'r', including start depth
r = c * time_vector / 2 + p.startdepth_m;

% Perform scan conversion

[cartesianImage, x_axis, y_axis] = scanconvert(harmonicIQ, r, theta);

logCompressedImage = 20 * log10(abs(cartesianImage));

% Define dynamic range and gain
dyn = 60; % Dynamic range in dB
gain = -90; 

% Display the log-compressed image
figure;
imagesc(logCompressedImage);
caxis([-dyn 0] - gain); % Adjust contrast and brightness
colormap(gray); % Set colormap to gray for ultrasound images
colorbar;

xlabel('Lateral Distance');
ylabel('Depth');
title('Log-Compressed Ultrasound Image Harmonic');



function demodulatedSignal = complexDemodulation(rfData, f0, fs, lowPassFilter)
    % rfData: The beamformed RF data matrix
    % f0: The transmit center frequency
    % fs: The sampling frequency
    % lowPassFilter: The low-pass filter coefficients or design

    % Time vector for the length of the RF data
    t = (0:size(rfData,1)-1)'/fs;

    % Downmixing
    % Create a complex exponential with frequency -f0 to shift the spectrum to baseband
    complexExponential = exp(-1i*2*pi*f0*t);
    downMixedSignal = rfData .* complexExponential;

    % Low-pass Filtering
    % Apply the low-pass filter to the downmixed signal
    demodulatedSignal = filtfilt(lowPassFilter, 1, downMixedSignal);
end

function demodulatedSignalHarmonic = complexDemodulationHarmonic(rfData, f0, fs)
    % rfData: The beamformed RF data matrix
    % f0: The transmit center frequency
    % fs: The sampling frequency

    % Generate time vector for the length of the RF data
    t = (0:size(rfData,1)-1)'/fs;

    % Downmixing to shift the spectrum centered around 2*f0 to baseband
    complexExponentialHarmonic = exp(-1i*2*pi*2*f0*t);
    downMixedSignal = rfData .* complexExponentialHarmonic;

    % Design a filter to suppress the fundamental frequency and pass the harmonic
    bpFilt = fir1(20, [1.8*f0 / (fs/2) 2.2*f0 / (fs/2)],"bandpass");

    % Apply the filter to the downmixed signal
    demodulatedSignalHarmonic = filtfilt(bpFilt, 1, downMixedSignal);
end



