% Load data
invivo_data = load("invivoData.mat");
load("simdata.mat");

% Task 1
% Constants
c = 1540; % Speed of sound in m/s
z = 0.01; % Depth in meters
x = 0; % Lateral position in meters

% Calculate Time of Flight (TOF) for simulated data
TOF_sim = TOF(z, c, x, elPosX);

% Calculate channel numbers and create x-axis
[sim_rows, sim_cols] = size(RFdata);
sim_xaxis = linspace(1, sim_cols, sim_cols);

% Visualize RF data on a logarithmic scale
figure;
imagesc('XData', sim_xaxis, 'YData', RF_t, 'CData', 20*log10(abs(RFdata)));
colormap('gray');
colorbar;
xlabel('Channel number');
title('RF data on logarithmic scale');
ylabel('Time (µs)');
hold on;
plot(TOF_sim, 'r', 'LineWidth', 2); % Plot TOF versus channel data
legend('TOF');

% Define coordinates
x = linspace(-2e-2, 2e-2, 256); % x-coordinates
z = linspace(0, 4.5e-2, 512).'; % z-coordinates
[X, Z] = meshgrid(x, z); % Create X and Z matrices

% Calculate TOF matrix
sim_elpos = permute(elPosX, [3 2 1]);
TOF_sim_mat = TOF(z, c, x, sim_elpos);

% Task 2

gain = -90;
dyn = 50;

% Interpolate TOF and visualize beamformed RF data
delayed_data = interpTOF(RFdata, RF_t, TOF_sim_mat);
bf_sim_rf = squeeze(sum(delayed_data, 1));
figure;
imagesc(20*log10(abs(bf_sim_rf)));
caxis([-dyn 0]-gain);
colormap("gray");
xlabel('Lateral Position'); % Adding X-axis label
ylabel('Depth'); % Adding Y-axis label
title('Beamformed RF Data on Logarithmic Scale'); % Adding title

% Calculate envelope and plot
sim_envelope = abs(hilbert(bf_sim_rf));
figure;
plot(sim_envelope);
xlabel('Sample Number'); % Adding X-axis label
ylabel('Envelope Amplitude'); % Adding Y-axis label
title('RF Data Envelope'); % Adding title

% Apodization for middle 32 elements
num_channels = size(delayed_data, 1);
middle_start = floor((num_channels - 32) / 2) + 1;
middle_end = middle_start + 32 - 1;
apod_32_middle = zeros(num_channels, 1);
apod_32_middle(middle_start:middle_end) = hamming(32);
apod_32_middle_reshaped = reshape(apod_32_middle, [length(apod_32_middle), 1, 1]);
apodized_data_middle32 = delayed_data .* apod_32_middle_reshaped;
bf_sim_rf_apod2 = squeeze(sum(apodized_data_middle32, 1));
figure;
imagesc(20*log10(abs(bf_sim_rf_apod2)));
colormap('gray');
caxis([-dyn 0]-gain);
xlabel('Lateral Position'); % Adding X-axis label
ylabel('Depth'); % Adding Y-axis label
title('Apodized Beamformed RF Data (Middle 32 Elements)'); % Adding title

% Apodization for first 32 elements
apod_32_first = zeros(size(delayed_data, 1), 1);
apod_32_first(1:32) = hamming(32);
apod_32_first_reshaped = reshape(apod_32_first, [length(apod_32_first), 1, 1]);
apodized_data_first32 = delayed_data .* apod_32_first_reshaped;
bf_sim_rf_apod1 = squeeze(sum(apodized_data_first32, 1));
dyn = 80;
gain = -60;
figure;
imagesc(20*log10(abs(bf_sim_rf_apod1)));
colormap('gray');
caxis([-dyn 0]-gain);
xlabel('Lateral Position'); % Adding X-axis label
ylabel('Depth'); % Adding Y-axis label
title('Apodized Beamformed RF Data (First 32 Elements)'); % Adding title



% Task 3: Apodization Visualization

% Define parameters
xpos = 0;
fnumber = 1;
z = linspace(0, 0.15, 512);

% Generate apodization matrix for F-number = 1
apod = generateApod(elPosX, xpos, z, fnumber);

% Visualize the apodization function versus depth and channel number
figure;
imagesc(1:length(elPosX), z, apod.');
xlabel('Channel number');
ylabel('Depth (m)');
title('Apodization Function vs. Depth and Channel Number for F-number = 1');
colormap('jet');
colorbar;

% Define a few key x-positions for demonstration
x_positions = [min(elPosX), 0, max(elPosX)]; % Start, middle, and end positions
f_numbers = [0.5, 1.5, 4]; % Low, suitable, and high F-number for comparison
z = linspace(0, 4.5e-2, 512); % Depth matrix for visualization

% Loop over selected F-numbers
for f_idx = 1:length(f_numbers)
    fnumber = f_numbers(f_idx);
    
    % Create a figure for each F-number
    figure;
    
    % Loop over selected x-positions
    for xpos_idx = 1:length(x_positions)
        xpos = x_positions(xpos_idx);
        
        % Generate the apodization matrix for this xpos and fnumber
        apod = generateApod(elPosX, xpos, z, fnumber);
        
        % Plot
        subplot(1, length(x_positions), xpos_idx);
        imagesc(1:length(elPosX), z*1e3, apod.'); % z in mm for better readability
        xlabel('Channel number');
        ylabel('Depth (mm)');
        title(sprintf('xpos=%.2f m, F#=%.1f', xpos, fnumber));
        colormap('jet');
        colorbar;
    end
    
    % Add a title for each F-number figure
    sgtitle(sprintf('Apodization Patterns for F-number = %.1f', fnumber));
end

% When the F-number increases, it becomes easier to maintain constant
% F-number with increasing distance. This is because the aperture is able to
% expand for longer.
% A large F-number means a narrower starting aperture that can focus deeper
% into the tissue. It increases the depth of field, but reduces the lateral
% resolution due to smaller number of elements contributing to the beam.







%Task 4
viv_elpos = invivo_data.elPosX;
viv_Rf = invivo_data.RFdata;
viv_Rf_t = invivo_data.RF_t;

x = linspace(-2e-2, 2e-2, 256); % x-coordinates
z = linspace(0, 4.5e-2, 512).'; % z-coordinates
[X, Z] = meshgrid(x, z); % Create X and Z matrices

viv_elpos = permute(viv_elpos, [3 2 1]);
TOF_viv = TOF(z, c, x, viv_elpos);

% Interpolate RF data based on TOF
delayed_data_viv = interpTOF(viv_Rf, viv_Rf_t, TOF_viv);

% Select an F-number for the expanding aperture
F = 1.2;
% Initialize expanding aperture matrix
expand_apod_viv = zeros(length(viv_elpos), length(z), length(x));

% Generate expanding aperture weights
for i = 1:length(x)
    expand_apod_viv(:,:,i) = generateApod(viv_elpos, x(i), z, F);
end

% Apply beamforming with expanding aperture
bf_viv_rf_ea = squeeze(sum(delayed_data_viv .* expand_apod_viv, 1));

% Visualization parameters
dyn = 70; % Dynamic range in dB
gain = -113; % Gain in dB

% Visualize the beamformed image
figure;
imagesc(x * 1e3, z * 1e3, 20 * log10(abs(bf_viv_rf_ea)));
colormap('gray');
caxis([-dyn, 0] - gain);
xlabel('Lateral position (mm)');
ylabel('Depth (mm)');
title('Ultrasound Image of Carotid Bifurcation');
axis equal;
colorbar;



% Define TOF function
function t = TOF(z, c, x, elposx)
    t_tx = z/c;
    t_rx = (sqrt((x-elposx).^2 + z.^2))/c;
    t = t_tx + t_rx;
end


