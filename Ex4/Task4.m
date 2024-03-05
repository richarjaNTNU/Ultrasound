%Load files
nearFieldPulse = load("NearFieldPulse.mat");
profile_AZ_7cm_full = load("profile_AZ_7cm_Full.mat");
profile_AZ_7cm_reduced = load("profile_AZ_7cm_reduced.mat");
profile_EL_7cm_full = load("profile_EL_7cm_Full.mat");
profile_EL_7cm_reduced = load("profile_EL_7cm_reduced.mat");
pulse_7cm_0dB_10x_aver = load("pulse_7cm_0dB_10x_aver.mat");
pulse_7cm_m8dB_10x_aver = load("pulse_7cm_-8dB_10x_aver.mat");
pulse_7cm_m14dB_10x_aver = load("pulse_7cm_-14dB_10x_aver.mat");
pulse_7cm_m14dB_no_aver = load("pulse_7cm_-14dB_no_aver.mat");
z_axis_12dB = load("z_axis_12dB.mat");

%Store RF files
pulse14_rf_noA = pulse_7cm_m14dB_no_aver.RF;
pulse14_rf_10A = pulse_7cm_m14dB_10x_aver.RF;
pulse8_rf_10A = pulse_7cm_m8dB_10x_aver.RF;
pulse0_rf_10A = pulse_7cm_0dB_10x_aver.RF;

%Time vectors   
time_pulse14_noA = linspace(pulse_7cm_m14dB_no_aver.Time,length(pulse14_rf_noA)*pulse_7cm_0dB_10x_aver.Ts,length(pulse14_rf_noA));
time_pulse14_10A = linspace(pulse_7cm_m14dB_10x_aver.Time,length(pulse14_rf_10A)*pulse_7cm_m14dB_10x_aver.Ts,length(pulse14_rf_10A));
time_pulse8_10A = linspace(pulse_7cm_m8dB_10x_aver.Time,length(pulse8_rf_10A)*pulse_7cm_m8dB_10x_aver.Ts,length(pulse8_rf_10A));
time_pulse0_10A = linspace(pulse_7cm_0dB_10x_aver.Time,length(pulse0_rf_10A)*pulse_7cm_0dB_10x_aver.Ts,length(pulse0_rf_10A));


%The lenght of the pulses are all the same.
num_samples = length(pulse14_rf_noA);
signal_stop_sig = idivide(num_samples,int16(4));
start_sig = 250;
noise1_start = 1;
noise1_end = 249;
noise2_start = signal_stop_sig + 1;
noise2_end = num_samples;

%Task 1 Pulse shapes
figure(2);
subplot(2,2,1),plot(time_pulse14_noA,pulse14_rf_noA),
title('Pulse 7cm -14dB no average'),
xlabel("\mus"),
ylabel("kPa"),

subplot(2,2,2), plot(time_pulse14_10A,pulse14_rf_10A),
title('Pulse 7cm -14dB 10x average'),
xlabel("\mus"),
ylabel("kPa"),

subplot(2,2,3),plot(time_pulse8_10A,pulse8_rf_10A),
title('Pulse 7cm -8dB 10x average'),
xlabel("\mus"),
ylabel("kPa"),

subplot(2,2,4),plot(time_pulse0_10A,pulse0_rf_10A),
title('Pulse 7cm 0dB 10x average'),
xlabel("\mus"),
ylabel("kPa"),


%Task 2 Fequency spectra, signal to noise ratio
%Signal spectra
figure(3), subplot(2,2,1),
plot(time_pulse14_noA(start_sig:signal_stop_sig),pulse14_rf_noA(start_sig:signal_stop_sig));
title('Pulse 7cm -14dB no average signal only');
xlabel("\mus");
ylabel("kPa");

subplot(2,2,2);
plot(time_pulse14_10A(start_sig:signal_stop_sig),pulse14_rf_10A(start_sig:signal_stop_sig));
title('Pulse 7cm -14dB 10x average signal only');
xlabel("\mus");
ylabel("kPa");

subplot(2,2,3);
plot(time_pulse8_10A(start_sig:signal_stop_sig),pulse8_rf_10A(start_sig:signal_stop_sig));
title('Pulse 7cm -8dB 10x average signal only');
xlabel("\mus");
ylabel("kPa");

subplot(2,2,4);
plot(time_pulse0_10A(start_sig:signal_stop_sig),pulse0_rf_10A(start_sig:signal_stop_sig));
title('Pulse 7cm 0dB 10x average signal only');
xlabel("\mus");
ylabel("kPa");

%Noise spectra
figure(4);
subplot(2,2,1);
plot(time_pulse14_noA([noise1_start:noise1_end noise2_start:noise2_end]),pulse14_rf_noA([noise1_start:noise1_end noise2_start:noise2_end]));
title('Pulse 7cm -14dB no average noise only');
xlabel("\mus");
ylabel("kPa");

subplot(2,2,2);
plot(time_pulse14_10A([noise1_start:noise1_end noise2_start:noise2_end]),pulse14_rf_10A([noise1_start:noise1_end noise2_start:noise2_end]));
title('Pulse 7cm -14dB 10x average noise only');
xlabel("\mus");
ylabel("kPa");

subplot(2,2,3);
plot(time_pulse8_10A([noise1_start:noise1_end noise2_start:noise2_end]),pulse8_rf_10A([noise1_start:noise1_end noise2_start:noise2_end]));
title('Pulse 7cm -8dB 10x average noise only');
xlabel("\mus");
ylabel("kPa");

subplot(2,2,4);
plot(time_pulse0_10A([noise1_start:noise1_end noise2_start:noise2_end]),pulse0_rf_10A([noise1_start:noise1_end noise2_start:noise2_end]));
title('Pulse 7cm 0dB 10x average noise only');
xlabel("\mus");
ylabel("kPa");


Nfft = 2048;
fs = 1/pulse_7cm_0dB_10x_aver.Ts;

%Signal spectra fft
freqtab = linspace(-0.5, 0.5, Nfft+1)*fs; freqtab(end) = [];
pow_spect_pulse14_rf_noA_sig = 20*log10( abs( fftshift( fft( pulse14_rf_noA(start_sig:signal_stop_sig), Nfft) ) ) );
pow_spect_pulse14_rf_10A_sig = 20*log10( abs( fftshift( fft( pulse14_rf_10A(start_sig:signal_stop_sig), Nfft) ) ) );
pow_spect_pulse8_rf_10A_sig = 20*log10( abs( fftshift( fft( pulse8_rf_10A(start_sig:signal_stop_sig), Nfft) ) ) );
pow_spect_pulse0_rf_10A_sig = 20*log10( abs( fftshift( fft( pulse0_rf_10A(start_sig:signal_stop_sig), Nfft) ) ) );

%Noise spectra fft
pow_spect_pulse14_rf_noA_noise  = 20*log10( abs( fftshift( fft( pulse14_rf_noA([noise1_start:noise1_end noise2_start:noise2_end]), Nfft) ) ) );
pow_spect_pulse14_rf_10A_noise = 20*log10( abs( fftshift( fft( pulse14_rf_10A([noise1_start:noise1_end noise2_start:noise2_end]), Nfft) ) ) );
pow_spect_pulse8_rf_10A_noise  = 20*log10( abs( fftshift( fft( pulse8_rf_10A([noise1_start:noise1_end noise2_start:noise2_end]), Nfft) ) ) );
pow_spect_pulse0_rf_10A_noise  = 20*log10( abs( fftshift( fft( pulse0_rf_10A([noise1_start:noise1_end noise2_start:noise2_end]), Nfft) ) ) );

% Limit the x-axis to 0-50 MHz
xlim_limit = 50e6;  % Set the x-axis limit to 50 MHz

%Plot spectra, blue for signal red for noise
figure(5);

subplot(2,2,1);
plot(freqtab, pow_spect_pulse14_rf_noA_sig, 'b', 'LineWidth', 0.5);
hold on;
plot(freqtab, pow_spect_pulse14_rf_noA_noise, 'r', 'LineWidth', 0.5);
title('Pulse 7cm -14dB Spectrum');
xlabel("Hz");
ylabel("dB");
legend('Signal', 'Noise');
xlim([0, xlim_limit]);  % Set x-axis limit
hold off;

subplot(2,2,2);
plot(freqtab, pow_spect_pulse14_rf_10A_sig, 'b', 'LineWidth', 0.5);
hold on;
plot(freqtab, pow_spect_pulse14_rf_10A_noise, 'r', 'LineWidth', 0.5);
title('Pulse 7cm -14dB 10x average Spectrum');
xlabel("Hz");
ylabel("dB");
legend('Signal', 'Noise');
xlim([0, xlim_limit]);  % Set x-axis limit
hold off;

subplot(2,2,3);
plot(freqtab, pow_spect_pulse8_rf_10A_sig, 'b', 'LineWidth', 0.5);
hold on;
plot(freqtab, pow_spect_pulse8_rf_10A_noise, 'r', 'LineWidth', 0.5);
title('Pulse 7cm -8dB 10x average Spectrum');
xlabel("Hz");
ylabel("dB");
legend('Signal', 'Noise');
xlim([0, xlim_limit]);  % Set x-axis limit
hold off;

subplot(2,2,4);
plot(freqtab, pow_spect_pulse0_rf_10A_sig, 'b', 'LineWidth', 0.5);
hold on;
plot(freqtab, pow_spect_pulse0_rf_10A_noise, 'r', 'LineWidth', 0.5);
title('Pulse 7cm 0dB 10x average Spectrum');
xlabel("Hz");
ylabel("dB");
legend('Signal', 'Noise');
xlim([0, xlim_limit]);  % Set x-axis limit


%Find max amplitude of spectra
[pow_spect_pulse14_rf_noA_sig_max,idx1] = max(pow_spect_pulse14_rf_noA_sig);
[pow_spect_pulse14_rf_10A_sig_max, idx2] = max(pow_spect_pulse14_rf_10A_sig);
[pow_spect_pulse8_rf_10A_sig_max,idx3] = max(pow_spect_pulse8_rf_10A_sig);
[pow_spect_pulse0_rf_10A_sig_max, idx4] = max(pow_spect_pulse0_rf_10A_sig);
fprintf('Max amplitude of spectrum -14dB no average is %4.2f at frequency %8.2f Hz',pow_spect_pulse14_rf_noA_sig_max, abs(freqtab(idx1)));
fprintf('Max amplitude of spectrum -14dB 10x average is %4.2f at frequency %8.2f Hz',pow_spect_pulse14_rf_10A_sig_max, abs(freqtab(idx2)));
fprintf('Max amplitude of spectrum -8dB 10x average is %4.2f at frequency %8.2f Hz',pow_spect_pulse8_rf_10A_sig_max, abs(freqtab(idx3)));
fprintf('Max amplitude of spectrum 0dB 10x average is %4.2f at frequency %8.2f Hz',pow_spect_pulse0_rf_10A_sig_max, abs(freqtab(idx4)));


%Task 3
profile_AZ_7cm_full_rf = profile_AZ_7cm_full.RF;
profile_AZ_7cm_full_pos = profile_AZ_7cm_full.Position;
profile_AZ_7cm_full_ts = profile_AZ_7cm_full.Ts;

profile_EL_7cm_full_rf = profile_EL_7cm_full.RF;
profile_EL_7cm_full_pos = profile_EL_7cm_full.Position;
profile_EL_7cm_full_ts = profile_EL_7cm_full.Ts;

AZ_intensity = sum(profile_AZ_7cm_full_rf.^2);
EL_intensity = sum(profile_EL_7cm_full_rf.^2);

AZ_intensity_db = 10*log10(AZ_intensity);
EL_intensity_db = 10*log10(EL_intensity);

max_int_AL_db = max(AZ_intensity_db);
max_int_EL_db = max(EL_intensity_db);

fprintf("Maximum intensity for AL %4.2f dB, and for EL %8.2f dB",max_int_AL_db,max_int_EL_db);

beam_AL_end_list = find(AZ_intensity_db < (max_int_AL_db-6));
beam_AL_end_idx = beam_AL_end_list(61);

beam_AL_end_list = find(AZ_intensity_db < (max_int_AL_db - 6));
beam_AL_end_idx = beam_AL_end_list(end);

beam_width_az_mm = abs(profile_AZ_7cm_full_pos(1, beam_AL_end_idx) - profile_AZ_7cm_full_pos(1, beam_AL_end_list(1)));

% Calculate two-sided beam width for elevation
beam_width_el_mm = abs(profile_EL_7cm_full_pos(2, beam_EL_end_list(end)) - profile_EL_7cm_full_pos(2, beam_EL_end_list(1)));

fprintf('Azimuth Beam Width: %.2f mm\n', beam_width_az_mm);
fprintf('Elevation Beam Width: %.2f mm\n', beam_width_el_mm);

figure(6);
subplot(2,1,1);
plot(profile_AZ_7cm_full_pos(1,:),AZ_intensity_db-max(AZ_intensity_db(:)));
title("AZ intensity with x position");
ylabel("intensity");
xlabel("mm");

subplot(2,1,2);
plot(profile_AZ_7cm_full_pos(2,:),AZ_intensity_db-max(AZ_intensity_db(:)));
title("AZ intensity with y position");
ylabel("intensity");
xlabel("mm");

figure(7);
subplot(2,1,1);
plot(profile_EL_7cm_full_pos(1,:),EL_intensity_db-max(EL_intensity_db(:)));
title("EL intensity with x position");
ylabel("intensity");
xlabel("mm");

subplot(2,1,2);
plot(profile_EL_7cm_full_pos(2,:),EL_intensity_db-max(EL_intensity_db(:)));
title("EL intensity with y position");
ylabel("intensity");
xlabel("mm");

%Task 6
load('profile_Az_7cm_full.mat');
timescale=Position(4,1)+Ts(1)*[1:size(RF,1)]; %Ts= Sampling rate
speedofsound=1500;%m/s
depthscale=timescale*speedofsound*1e3;%mm
azimutpositionscale=Position(1,:)*1e3; %mm
imagesc(azimutpositionscale,depthscale,RF); %Neg. pressure dark
colormap(gray);

subplot(2,1,1);
load('profile_AZ_7cm_reduced.mat');
timescale=Position(4,1)+Ts(1)*[1:size(RF,1)]; %Ts= Sampling rate
speedofsound=1500;%m/s
depthscale=timescale*speedofsound*1e3;%mm
azimutpositionscale=Position(1,:)*1e3; %mm
imagesc(azimutpositionscale,depthscale,RF); %Neg. pressure dark
colormap(gray); 

subplot(2,1,2);
load('profile_Az_7cm_full.mat');
timescale=Position(4,1)+Ts(1)*[1:size(RF,1)]; %Ts= Sampling rate
speedofsound=1500;%m/s
depthscale=timescale*speedofsound*1e3;%mm
azimutpositionscale=Position(1,:)*1e3; %mm
imagesc(azimutpositionscale,depthscale,RF); %Neg. pressure dark
colormap(gray); 

