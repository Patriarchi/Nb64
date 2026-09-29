%Created by Paul J. Lamothe-Molina
%Version 1.1 - 20250213
%Script developed to analyse photometry recordings of nLightG2 in dorsal HPC after GqDREADD activation of LC. 
%Script developed for doric photometry files. 

close all
clear all

%% Plot colors and definitions:
% Number of colors in the colormap
numColors = 256;
pink = [230, 114, 191] / 255; % RGB for pink (divided by 255 to normalize to [0, 1])
white = [255, 255, 255] / 255; % RGB for white (divided by 255 to normalize to [0, 1])
green = [0, 127, 0] / 255; % RGB for green (divided by 255 to normalize to [0, 1])
plum = [212, 161, 202] / 255; % RGB for plum (divided by 255 to normalize to [0, 1])

% Linearly interpolate between white and pink
for i = 1:numColors
    WhiteToPink(i, :) = (i-1) / (numColors-1) * pink + (1 - (i-1) / (numColors-1)) * white;
end

% Linearly interpolate between white and green
for i = 1:numColors
    WhiteToGreen(i, :) = (i-1) / (numColors-1) * green + (1 - (i-1) / (numColors-1)) * white;
end



%% Signal extraction from Doric files

%Step 1: set path to code and select directory where the .doric file is
addpath ('K:\PLM\Code\Matlab\') % add path with the matalb functions
addpath ('K:\PLM\Code\Matlab\Toolboxes\Matlab Functions_MH_ZK') % add path with the matalb functions

%give the folder with the raw data
DataFolder = uigetdir('D:\Analysis\ANNI\in_vivo\GqDREADD_LC_dHPC\mouseLR_M_ac\CFC'); %give the folder with the raw data
cd (DataFolder); % set the folder for saving the data
list_of_files = dir ('*.doric'); %reads name and other characteristics of all the .doric files in that directory
filename = list_of_files.name;  %defines the filename variable based on the given file name

%Step 2: Get information about the Doric file:
DoricInfo = h5info(filename);
%read data recorded from defferent channels
%create cells with datasets
%first colum: name of the channel
for Event_Number = 1:size (DoricInfo.Groups.Groups.Groups,1)
    DoricData{Event_Number,1}=DoricInfo.Groups.Groups.Groups(Event_Number).Datasets.Name ;
    CurrentChannel=sprintf('%s%s%s', DoricInfo.Groups.Groups.Groups(Event_Number).Name,'/',DoricInfo.Groups.Groups.Groups(Event_Number).Datasets.Name );
    DoricData{Event_Number,2}=h5read(filename,CurrentChannel);
end

InjectionTime = 180; %Injection time in seconds

%Step 3: Extract the signals:
Isoemissive_405nm_Demodulated = DoricData{1,2};
nLightG2_Demodulated = DoricData{2,2};
SWD_RawSignal = DoricData{3,2};
Shock_TTL = DoricData{6,2}; %Shock TTL channel
%recording to the 
ConsoleTime = DoricData{7,2};

%% Signal processing

%Step 1: clean signal from NaN
nLightG2_Signal_Demodulated_Clean = nLightG2_Demodulated;
Isoemissive_405nm_Demodulated_Clean = Isoemissive_405nm_Demodulated;
Shock_TTL_Clean = Shock_TTL;
ConsoleTime_Clean = ConsoleTime;
SWD_RawSignal_Clean = SWD_RawSignal;

k=1; % loop that detects NaN values of both nLightG2 and Isoemissive signals and created a matrix of thoese indices. 
for Event_Number=1:length(Isoemissive_405nm_Demodulated)
if isnan(nLightG2_Demodulated(Event_Number))==1||isnan(Isoemissive_405nm_Demodulated(Event_Number))==1
   
    NaNIndex2Delete(k) = Event_Number;
    k = k+1;
end
end


%Removes those NaN indicies from all channels
nLightG2_Signal_Demodulated_Clean(NaNIndex2Delete) = [];
Isoemissive_405nm_Demodulated_Clean(NaNIndex2Delete) = [];
ConsoleTime_Clean(NaNIndex2Delete) = [];
Shock_TTL_Clean(NaNIndex2Delete) = [];
SWD_RawSignal_Clean(NaNIndex2Delete) = [];


%Step 2: plot the data of the whole session
Y_Shifted_405 = Isoemissive_405nm_Demodulated_Clean + (mean(nLightG2_Signal_Demodulated_Clean));
Y_Shifted_TTL = Shock_TTL_Clean + (mean(nLightG2_Signal_Demodulated_Clean));

NamePlot='Raw nLightG2b Signal and Isoemissive';
figure;
hold on;
plot(ConsoleTime_Clean, nLightG2_Signal_Demodulated_Clean, 'Color', patlab_colors('DarkGreen'), 'DisplayName', 'nLightG2b Raw');
plot(ConsoleTime_Clean, Isoemissive_405nm_Demodulated_Clean, 'Color', patlab_colors('Plum'), 'DisplayName', 'Isoemissive');
plot(ConsoleTime_Clean, Y_Shifted_TTL, 'Color', patlab_colors('FireBrick'), 'DisplayName', 'Shock TTL');
title(['Raw nLightG2 and Isoemissive of ', filename], 'Interpreter', 'none');
xlabel 'Time (s)';
ylabel 'Detector voltage (V)';
legend('show');
hold off;
SavePlots(NamePlot);



%% Downsampling by rolling mean - Only run for 12000 samples/s files
dT_Doric_average=max(ConsoleTime)/length(nLightG2_Signal_Demodulated_Clean);
FpS = 1/dT_Doric_average; %Frames per Second
DownSample = 100;

% Length of the downsampled signal
new_length = floor(length(ConsoleTime_Clean) / DownSample);
% New sampling rate after downsampling
FpS_downsampled = FpS / DownSample;
% Create the new time vector based on the new sampling rate
ConsoleTime_Clean_Downsampled = (0:new_length-1) / FpS_downsampled;

% Define the cutoff frequency (e.g., half of the new sample rate)
cutoff_freq = FpS_downsampled / 2;  % Nyquist frequency of the downsampled signal
nyquist_freq = FpS / 2;  % Nyquist frequency of the original signal

% Design a low-pass filter (Butterworth filter as an example)
[b, a] = butter(4, cutoff_freq / nyquist_freq);  % 4th-order Butterworth filter

% Apply the filter to the original signal
filtered_signal_nLightG2 = filter(b, a, nLightG2_Signal_Demodulated_Clean);
filtered_signal_405nm = filter(b, a, Isoemissive_405nm_Demodulated_Clean);


% Initialize the averaged signal
Isoemissive_405nm_Demodulated_Clean_Downsampled = zeros(1, floor(length(filtered_signal_405nm)/DownSample));

% Loop to average every 100 samples
for i = 1:length(Isoemissive_405nm_Demodulated_Clean_Downsampled)
    start_idx = (i-1) * DownSample + 1;
    end_idx = i * DownSample;
    Isoemissive_405nm_Demodulated_Clean_Downsampled(i) = mean(filtered_signal_405nm(start_idx:end_idx));
end

% Initialize the averaged signal
nLightG2_Signal_Demodulated_Clean_Downsampled = zeros(1, floor(length(filtered_signal_nLightG2)/DownSample));

% Loop to average every 100 samples
for i = 1:length(nLightG2_Signal_Demodulated_Clean_Downsampled)
    start_idx = (i-1) * DownSample + 1;
    end_idx = i * DownSample;
    nLightG2_Signal_Demodulated_Clean_Downsampled(i) = mean(filtered_signal_nLightG2(start_idx:end_idx));
end

nLightG2_Signal_Demodulated_Clean_Downsampled = nLightG2_Signal_Demodulated_Clean_Downsampled';
Isoemissive_405nm_Demodulated_Clean_Downsampled = Isoemissive_405nm_Demodulated_Clean_Downsampled';
ConsoleTime_Clean_Downsampled = ConsoleTime_Clean_Downsampled';

NamePlot='Downsampled Raw nLightG2b Signal and Isoemissive';
figure;
hold on;
plot(ConsoleTime_Clean_Downsampled, nLightG2_Signal_Demodulated_Clean_Downsampled, 'Color', green, 'DisplayName', 'nLightG2b Raw');
plot(ConsoleTime_Clean_Downsampled, Isoemissive_405nm_Demodulated_Clean_Downsampled, 'Color', plum, 'DisplayName', 'Isoemissive');
title(['Downsampled Raw nLightG2 and Isoemissive of ', filename], 'Interpreter', 'none');
%plot(ConsoleTime_Clean, Y_Shifted_TTL, 'c', 'DisplayName', '488 nm laser TTL');
xlabel 'Time (s)';
ylabel 'Detector voltage (V)';
legend('show');
hold off;
SavePlots(NamePlot);

% %% Spectral Analysis
% 
% % Define parameters
% L = length(SWD_RawSignal_Clean); % Length of the signal
% t = (0:L-1) / FpS; % Time vector
% 
% % Compute the Fourier Transform
% Y = fft(SWD_RawSignal_Clean);
% P2 = abs(Y / L); % Two-sided spectrum
% P1 = P2(1:floor(L/2)+1); % Single-sided spectrum
% P1(2:end-1) = 2 * P1(2:end-1); % Normalize amplitude
% 
% % Frequency vector
% f = FpS * (0:(L/2)) / L;
% 
% % Plot power spectrum
% figure;
% plot(f, P1, 'b', 'LineWidth', 1.5);
% xlabel('Frequency (Hz)');
% ylabel('Power');
% title('Power Spectrum of SWD_RawSignal');
% grid on;
% xlim([0 FpS/2]); % Show only the positive half of frequencies
% 



%% Run this for low-sampled recordings
dT_Doric_average=max(ConsoleTime)/length(nLightG2_Signal_Demodulated_Clean);
FpS = 1/dT_Doric_average; %Frames per Second
ConsoleTime_Clean_Downsampled = ConsoleTime_Clean;
Isoemissive_405nm_Demodulated_Clean_Downsampled = Isoemissive_405nm_Demodulated_Clean;
nLightG2_Signal_Demodulated_Clean_Downsampled = nLightG2_Signal_Demodulated_Clean;
FpS_downsampled = FpS;
%%

FiveSecondIndex = 5*FpS_downsampled;
ConsoleTime_Clean_Downsampled = ConsoleTime_Clean_Downsampled(180*FpS_downsampled:end);
Isoemissive_405nm_Demodulated_Clean_Downsampled = Isoemissive_405nm_Demodulated_Clean_Downsampled(180*FpS_downsampled:end);
nLightG2_Signal_Demodulated_Clean_Downsampled = nLightG2_Signal_Demodulated_Clean_Downsampled(180*FpS_downsampled:end);
Shock_TTL_Clean_Truncated = Shock_TTL_Clean(180*FpS_downsampled:end);

%Step 3: signal filtering

%low-pass
f_lowpass = 5;
Isoemissive_405nm_Demodulated_Clean_Filtered=lowpass(Isoemissive_405nm_Demodulated_Clean_Downsampled, f_lowpass, FpS_downsampled);
nLightG2_Signal_Demodulated_Clean_Filtered=lowpass(nLightG2_Signal_Demodulated_Clean_Downsampled, f_lowpass, FpS_downsampled);


%remove edge artifact
FilteredSignal_Length = length(nLightG2_Signal_Demodulated_Clean_Filtered);
ConsoleTime_Clean_noEdge = ConsoleTime_Clean_Downsampled(100:(FilteredSignal_Length - 100), :);
Isoemissive_405nm_Demodulated_Clean_filtered_noEdge = Isoemissive_405nm_Demodulated_Clean_Filtered(100:(FilteredSignal_Length - 100), :);
nLightG2_Signal_Demodulated_Clean_Clean_Filtered_noEdge = nLightG2_Signal_Demodulated_Clean_Filtered(100:(FilteredSignal_Length - 100), :);
Shock_TTL_Clean_NoEdge = Shock_TTL_Clean_Truncated(100:(FilteredSignal_Length - 100), :);


%plot to with edge artifact
NamePlot = 'nLightG2 vs. Isoemissive RAW power - Edge Artifact';
figure;
plot(nLightG2_Signal_Demodulated_Clean_Filtered, Isoemissive_405nm_Demodulated_Clean_Filtered, 'o', 'Color', patlab_colors('OrangeRed'), 'DisplayName', 'With Artifact');
title(['nLightG2 vs. Isoemissive RAW power - Edge Artifact of ', filename], 'Interpreter', 'none', 'FontSize', 8);

xlabel('Isosbestic');
ylabel('Functional');
legend('show');
grid on;
SavePlots(NamePlot);

%plot to confirm removal of edges
NamePlot = 'nLightG2 vs. Isoemissive RAW power - No Edge Artifact';
figure;
plot(nLightG2_Signal_Demodulated_Clean_Clean_Filtered_noEdge, Isoemissive_405nm_Demodulated_Clean_filtered_noEdge, 'o', 'Color', patlab_colors('Turquoise'), 'DisplayName', 'Without Artifact');
title(['nLightG2 vs. Isoemissive RAW power - No Edge of ', filename], 'Interpreter', 'none', 'FontSize', 8);
xlabel('Isosbestic');
ylabel('Functional');
legend('show');
grid on;
SavePlots(NamePlot);

%% Correction to Isoemissive - Use when no photobleaching occurs or with traces without a baseline
dF_Div = (nLightG2_Signal_Demodulated_Clean_Clean_Filtered_noEdge ./ Isoemissive_405nm_Demodulated_Clean_filtered_noEdge)+1;
dF_Sub = nLightG2_Signal_Demodulated_Clean_Clean_Filtered_noEdge - Isoemissive_405nm_Demodulated_Clean_filtered_noEdge;

%plot for dF (DIV)
NamePlot = 'Functional Signal - Corrected by Isoemissive (DIV)';
figure;
plot(ConsoleTime_Clean_noEdge, dF_Div, 'Color', patlab_colors('DarkGreen'));
title(['nLightG2 dF (DIV) of ', filename], 'Interpreter', 'none', 'FontSize', 8);
xlabel('Time (s)');
ylabel('arbitrary units');
SavePlots(NamePlot);


%plot for dF (SUB)
NamePlot = 'Functional Signal - Corrected by Isoemissive (DIV)';
figure;
hold on;
plot(ConsoleTime_Clean_noEdge, dF_Sub, 'Color', patlab_colors('DarkGreen'));
plot(ConsoleTime_Clean_noEdge, Shock_TTL_Clean_NoEdge, 'Color', patlab_colors('FireBrick'), 'DisplayName', 'Shock TTL');
title(['nLightG2 dF (SUB) of ', filename], 'Interpreter', 'none', 'FontSize', 8);
xlabel('Time (s)');
ylabel('arbitrary units');
SavePlots(NamePlot);

%% PB correction polynomial
Y_405nm = Isoemissive_405nm_Demodulated_Clean_filtered_noEdge; % (5*FpS_downsampled:end);
Y_465nm = nLightG2_Signal_Demodulated_Clean_Clean_Filtered_noEdge; %(5*FpS_downsampled:end);
ConsoleTime_Truncated = ConsoleTime_Clean_noEdge; %(5*FpS_downsampled:end);

coefficients_LinFit_405nm = polyfit(Y_405nm, ConsoleTime_Truncated, 1);
Slope_405nm = coefficients_LinFit_405nm(1);
Intercept_405nm = coefficients_LinFit_405nm(2);
correction_vals_405nm = Isoemissive_405nm_Demodulated_Clean_filtered_noEdge*Slope_405nm + Intercept_405nm;
Isoemissive_405nm_PBcorrected_LinFit = Isoemissive_405nm_Demodulated_Clean_filtered_noEdge - correction_vals_405nm;

coefficients_LinFit_465nm = polyfit(Y_465nm, ConsoleTime_Truncated, 1);
Slope_465nm = coefficients_LinFit_465nm(1);
Intercept_465nm = coefficients_LinFit_465nm(2);
correction_vals_465nm = nLightG2_Signal_Demodulated_Clean_Clean_Filtered_noEdge*Slope_465nm + Intercept_465nm;
Functional_465nm_PBcorrected_LinFit = nLightG2_Signal_Demodulated_Clean_Clean_Filtered_noEdge - correction_vals_465nm;

figure;
plot(ConsoleTime_Clean_noEdge, correction_vals_465nm);

figure;
plot(ConsoleTime_Clean_noEdge, Functional_465nm_PBcorrected_LinFit);

%%
%simple linear regression between Isosbestic and Functional
X_405nm = Isoemissive_405nm_PBcorrected_LinFit;
Y_465nm = Functional_465nm_PBcorrected_LinFit;

coefficients_MotionCorrection_465nm = polyfit(X_405nm, Y_465nm, 1);
MotionCorrection_Slope_465nm = coefficients_MotionCorrection_465nm(1);
MotionCorrection_Intercept_465 = coefficients_MotionCorrection_465nm(2);
EstimatedMotionArtifact_465nm = MotionCorrection_Slope_465nm*X_405nm + MotionCorrection_Intercept_465;
nLightG2_PBcorrected_SUB_MOcorrected = (Functional_465nm_PBcorrected_LinFit - EstimatedMotionArtifact_465nm);



%change the size of the TTL to fit the plots.
%Y_10x_Truncated_TTL = Laser488nm_TTL_Clean_noEdge / 10;
%Y_100x_Truncated_TTL = Laser488nm_TTL_Clean_noEdge / 100;

%plot corrected signals with TTL
NamePlot='Motion and photobleaching corrected nLightG2 Signal and Laser TTL';
figure;
title('Motion and photobleaching corrected nLightG2 Signal and Laser TTL');
hold on;
plot(ConsoleTime_Clean_noEdge, nLightG2_PBcorrected_SUB_MOcorrected, 'r', 'DisplayName', 'nLightG2 Corrected');
xlabel 'Time (s)';
ylabel 'Detector voltage (V)';
legend('show');
hold off;
SavePlots(NamePlot);

%% dFF
% 
% Assuming nLightG2 is your input signal

% Compute the 10th percentile threshold
percentile5_threshold = prctile(dF_Sub, 5);

% Find values below or equal to the 10th percentile threshold
low_values = dF_Sub(dF_Sub <= percentile5_threshold);

% Calculate the baseline as the mean of those values
Baseline10PER_nLightG2 = mean(low_values);

% Compute dFF using the formula: (Fn - Baseline10PER_nLightG2) / Baseline10PER_nLightG2
%dFF_nLightG2 = 100*(dF_Sub - Baseline10PER_nLightG2) / Baseline10PER_nLightG2;

dF_Baseline = mean(dF_Sub(200:450*FpS_downsampled));
dFF_nLightG2 = 100*(dF_Sub - dF_Baseline) ./ dF_Baseline;
dFF_nLightG2_Smooth = medfilt1(dFF_nLightG2, 100);

dFF_nLightG2_5per = 100*(dF_Sub - Baseline10PER_nLightG2) ./ Baseline10PER_nLightG2;


NamePlot='dFF nLightG2_Signal';
figure;
title('dFF nLightG2_Signal by percentile');
hold on;
subplot(2, 1, 1);
hold on;
plot(ConsoleTime_Clean_noEdge, dFF_nLightG2, 'Color', patlab_colors('DarkGreen'), 'DisplayName', 'Functional Signal');
plot(ConsoleTime_Clean_noEdge, dFF_nLightG2_Smooth, 'Color', patlab_colors('PaleGreen'), 'DisplayName', 'Smooth');
plot(ConsoleTime_Clean_noEdge, Shock_TTL_Clean_NoEdge, 'Color', patlab_colors('FireBrick'), 'DisplayName', 'Shock TTL');
hold off;
xlabel 'Time (s)';
ylabel 'dFF %';
legend('show');
subplot(2, 1, 2);
plot(ConsoleTime_Clean_noEdge, dFF_nLightG2_5per, 'k', 'DisplayName', 'dFF (5 low percentile BL)');
xlabel 'Time (s)';
ylabel 'dFF %';
legend('show');
hold off;
SavePlots(NamePlot);



%% Z-Score calculation

%Z-Score Whole Trace

Mean4zScore_WholeTrace = mean(dF_Sub);
StDev4zScore_WholeTrace = std(dF_Sub);
zScore_nLightG2_WholeTrace = (dF_Sub - Mean4zScore_WholeTrace) ./ StDev4zScore_WholeTrace;
zScore_nLightG2_WholeTrace_Smooth = medfilt1(zScore_nLightG2_WholeTrace, 150);

Mean4zScore_Baseline = mean(dF_Sub(250*FpS:550*FpS_downsampled));
StDev4zScore_Baseline = std(dF_Sub(250*FpS:550*FpS_downsampled));
zScore_nLightG2_Baseline = (dF_Sub - Mean4zScore_Baseline) ./ StDev4zScore_Baseline;
zScore_nLightG2_Baseline_Smooth = medfilt1(zScore_nLightG2_Baseline, 150);



NamePlot='zScore nLightG2_Signal';
figure;
title(['zScore of nLightG2 of ', filename], 'Interpreter', 'none', 'FontSize', 8);
hold on;
subplot(2,1,1);
hold on;
plot(ConsoleTime_Clean_noEdge, zScore_nLightG2_WholeTrace, 'Color', patlab_colors('DarkGreen'), 'DisplayName', 'zScore Whole Trace');
plot(ConsoleTime_Clean_noEdge, zScore_nLightG2_WholeTrace_Smooth, 'Color', patlab_colors('LightGray'), 'DisplayName', 'zScore From Whole Trace smooth');
hold off;
xlabel 'Time (s)';
ylabel 'zScore';
legend('show');
subplot(2,1,2);
hold on;
plot(ConsoleTime_Clean_noEdge, zScore_nLightG2_Baseline, 'Color', patlab_colors('PaleGreen'), 'DisplayName', 'zScore From Baseline');
plot(ConsoleTime_Clean_noEdge, zScore_nLightG2_Baseline_Smooth, 'Color', patlab_colors('DarkGray'), 'DisplayName', 'zScore From Baseline smooth');
hold off;
xlabel 'Time (s)';
ylabel 'zScore';
legend('show');
hold off;
SavePlots(NamePlot);


%% Heatmaps
timeIndices = 0:600*FpS:3600*FpS;
xticks_values = 0:30:3600;

NamePlot='Z-Score Heatmap nLightG2';
figure;
imagesc(zScore_nLightG2_Baseline_Smooth');
title('Z-Score (Excluding Stim) Heatmap nLightG2')
colorbar;
xlabel('Time (s)');
xticks(timeIndices);
xticklabels(xticks_values);
colormap (WhiteToGreen);
%caxis([-2, 15]);
SavePlots(NamePlot);

%% Individual trials anaylsis

event_length_seconds = 2; %Footshock in seconds
num_trials = 10;
ITI = 125; % intertrial interval in seconds
Time2Sync = 600:ITI:600+((num_trials - 1)*ITI); %timestamps for shocks vector
event_timestamps = Time2Sync;
window_before = 20;
window_after = 60;
min_length = Inf;
TrialLength = window_before + window_after + event_length_seconds;
TimeInterval = 1/FpS;

ConsoleTime_Trials = -0:TimeInterval:TrialLength;
ConsoleTime_Trials = ConsoleTime_Trials - 20;
%%
%Analysis of individual events for z-score nLightG2

% Divide the dataset into subsections around each event
for i = 1:num_trials
    % Find the index corresponding to the event time
    [~, event_index] = min(abs(ConsoleTime_Clean_noEdge - event_timestamps(i)));

    % Calculate start and end indices based on event length and windows
    start_index = find(ConsoleTime_Clean_noEdge >= (event_timestamps(i) - window_before), 1);
    end_index = find(ConsoleTime_Clean_noEdge >= (event_timestamps(i) + event_length_seconds + window_after), 1);

    % Extract data for the current subsection
    data_subsections_nLightG2{i} = zScore_nLightG2_Baseline_Smooth(start_index:end_index, :);

        % Update min_length if the current subsection is shorter
    min_length = min(min_length, size(data_subsections_nLightG2{i}, 1));
end

% Initialize a matrix to store all subsection data
individual_events_nLightG2 = NaN(min_length, num_trials);

% Concatenate data from all subsections into a single matrix
for i = 1:num_trials
    individual_events_nLightG2(:, i) = data_subsections_nLightG2{i}(1:min_length, 1); % Only taking signal values
end

%%
%%
%Analysis of individual events for dF nLightG2

% Divide the dataset into subsections around each event
for i = 1:num_trials
    % Find the index corresponding to the event time
    [~, event_index] = min(abs(ConsoleTime_Clean_noEdge - event_timestamps(i)));

    % Calculate start and end indices based on event length and windows
    start_index = find(ConsoleTime_Clean_noEdge >= (event_timestamps(i) - window_before), 1);
    end_index = find(ConsoleTime_Clean_noEdge >= (event_timestamps(i) + event_length_seconds + window_after), 1);

    % Extract data for the current subsection
    data_subsections_dF_nLightG2{i} = dF_Sub(start_index:end_index, :);

        % Update min_length if the current subsection is shorter
    min_length = min(min_length, size(data_subsections_dF_nLightG2{i}, 1));
end

% Initialize a matrix to store all subsection data
individual_events_dF_nLightG2 = NaN(min_length, num_trials);

% Concatenate data from all subsections into a single matrix
for i = 1:num_trials
    individual_events_dF_nLightG2(:, i) = data_subsections_nLightG2{i}(1:min_length, 1); % Only taking signal values
end



%%
ConsoleTime_inMinutes = ConsoleTime_Clean_noEdge / 60;
%%

AUC_zScore_nLightG2 = trapz(abs(zScore_nLightG2_Baseline));

%% Save Variables
save([filename '.mat'])
