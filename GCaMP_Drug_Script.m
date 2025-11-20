%% Visualize Raw Data
raw_entry = rawTDT1("D:\Will\WD-035-FL-251117-115346\NA-FL15-251118-163915", 3, 0, 1, 'PC0 PC1');
%% Add to Raw Data Structure
raw_structure = packageFP3('ratopioidraw.mat', raw_entry, 1); % add individual entry with identifying info to a full structure
prompt = 'File Save Name: '; % sanity check to ensure the right file name gets saves
save_name = input(prompt);
save(save_name, "raw_structure", '-v7.3');
%% Visualize Raw Data
idx = 34; % index raw structure to visualize individul entry
row = raw_structure(idx);
time = row.Time_Vector; % extract trace data from raw structure
sensor = row.Sensor_Signal;
control = row.Control_Signal;
ttl_matrix = row.TTLs; % extract TTLs from raw structure
pc1 = ttl_matrix{2};
figure; % plot data
hold on;
line_sensor = plot(time, sensor, 'g', 'DisplayName', 'Sensor');
line_control = plot(time, control, 'm', 'DisplayName', 'Control');
title('Raw Data');
ylabel('Signal (mV)');
xlabel('Time (s)');
xline(pc1(1), '-', 'Injection'); % plots first TTL as the injection
legend([line_sensor line_control]);
%% Normalize GCaMP Data
idx = 50; % index raw structure for individual entries
row = raw_structure(idx);
ttl_matrix = row.TTLs; % extract TTL data
pc1 = ttl_matrix{2};
beginning = pc1(1); % set beginning of drug period
row.Time_Vector = smooth(row.Time_Vector, 0.002, 'lowess'); % smooth raw data with a local regression using a weighted linear least squares
row.Sensor_Signal = smooth(row.Sensor_Signal, 0.002, 'lowess');
row.Control_Signal = smooth(row.Control_Signal, 0.002, 'lowess');
row_downsampled = downsampleFP1(row, 25, 0); % use a custom function to downsample data to 25 Hz
time = row.Time_Vector;
ending = time(end); % set the timepoint for the end of the recording
baseline_range = [beginning-600 beginning]; % set the time range for the baseline
normalization_range = [beginning-600 ending]; % set the time range for data to be normalized
normalized = normalizeFP3(row_downsampled, 0, 'zscore', 1, 1, baseline_range, normalization_range); % normalize data to isosbestic, save and plot signal as z-score of ΔF/F 
%% Manually Adjust Baseline to 0
t = normalized.Time_Vector; % extract signal data from normalized recording
normalized_signal = normalized.Signal_Vector;
ttl_matrix = normalized.TTLs; % extract TTL data
pc1 = ttl_matrix{2};
injection = pc1(1); % find first TTL that represents injection
baseline_idx = find(t < injection); % index times before the injection to find the baseline
baseline = normalized_signal(baseline_idx);
baseline_time = t(baseline_idx);
mean_baseline = mean(baseline); % find the average of the baseline signal
figure; % plot the normalized signal before the adjusting the baseline
hold on;
plot(t, normalized_signal, 'b');
yline(mean_baseline);
yline(0); % show x-axis (y = 0)
xline(injection); % show where injection occurs
area(baseline_time, normalized_signal(baseline_idx), 0, 'FaceColor', 'r'); % shades the area between the baseline signal and 0
difference = 0 - mean_baseline; % finds difference between baseline average and 0
normalized_signal = normalized_signal + difference; % adds this difference to each data point to bring the average of the baseline period to 0
figure; % plot the new signal with the adjusted baseline
hold on;
plot(t, normalized_signal);
yline(0);
xline(injection);
area(baseline_time, normalized_signal(baseline_idx), 0, 'FaceColor', 'r');
normalized_adjusted = normalized; % create a seperate structure for this data to preserve the original
normalized_adjusted.Signal_Vector = normalized_signal; % adds the new adjusted signal to this new structure
%% Save Structure of Normalized Data
normalized_structure = packageFP3('gcampopioidnormfinal.mat', normalized_adjusted, 1); % add the adjusted structure, with identifiers, to the full structure of normalized data
prompt = 'File Save Name: '; % sanity check to ensure that the correct file name is used
save_name = input(prompt);
save(save_name, "normalized_structure");
%% AUC for GCaMPs
idx = 50;
row = normalized_structure(idx); % index the normalized data structure for individual entries
bucketed_auc = phaucFP1(row, 300, 3600, 'PC1', 1); % use custom function to find area under the curve in specified time buckets
%% Save Structure of AUC Data
auc_bucketed = packageFP3('gcampopioidaucfinal.mat', bucketed_auc, 1); % save AUC data to its own structure
prompt = 'File Save Name: ';
save_name = input(prompt);
save(save_name, "auc_bucketed");
%% Combine/Average Recordings
combined = combineFP3(auc_bucketed, [{'Location'},{'NA'}; {'Sensor'}, {'Ca'}; {'Experiment'}, {'Buprenorphine 1'}]); % use custom function to combine trace data for structure elements with specified identifiers
%% Save Structure of Combined Traces
combined_structure = packageFP3('gcampopioidcombined.mat', combined, 1); % add data to a structure of combined recordings
prompt = 'File Save Name: ';
save_name = input(prompt);
save(save_name, "combined_structure");
%% Plot Averaged Traces Together
saline = combined_structure(1); % combine pre and post injection average traces and SEM for each relevant condition
saline_pre_trace = saline.Signal_Pre_Timepoint_Averaged;
saline_pre_sem = saline.Signal_Pre_Timepoint_SEM;
saline_post_trace = saline.Signal_Post_Timepoint_Averaged;
saline_post_sem = saline.Signal_Post_Timepoint_SEM;
saline_combined = [saline_pre_trace saline_post_trace];
saline_sem = [saline_pre_sem saline_post_sem];
fnz003 = combined_structure(2);
fnz003_pre_trace = fnz003.Signal_Pre_Timepoint_Averaged;
fnz003_pre_sem = fnz003.Signal_Pre_Timepoint_SEM;
fnz003_post_trace = fnz003.Signal_Post_Timepoint_Averaged;
fnz003_post_sem = fnz003.Signal_Post_Timepoint_SEM;
fnz003_combined = [fnz003_pre_trace fnz003_post_trace];
fnz003_sem = [fnz003_pre_sem fnz003_post_sem];
fnz01 = combined_structure(3);
fnz01_pre_trace = fnz01.Signal_Pre_Timepoint_Averaged;
fnz01_pre_sem = fnz01.Signal_Pre_Timepoint_SEM;
fnz01_post_trace = fnz01.Signal_Post_Timepoint_Averaged;
fnz01_post_sem = fnz01.Signal_Post_Timepoint_SEM;
fnz01_combined = [fnz01_pre_trace fnz01_post_trace];
fnz01_sem = [fnz01_pre_sem fnz01_post_sem];
dfnz03 = combined_structure(4);
dfnz03_pre_trace = dfnz03.Signal_Pre_Timepoint_Averaged;
dfnz03_pre_sem = dfnz03.Signal_Pre_Timepoint_SEM;
dfnz03_post_trace = dfnz03.Signal_Post_Timepoint_Averaged;
dfnz03_post_sem = dfnz03.Signal_Post_Timepoint_SEM;
dfnz03_combined = [dfnz03_pre_trace dfnz03_post_trace];
dfnz03_sem = [dfnz03_pre_sem dfnz03_post_sem];
dfnz1 = combined_structure(5);
dfnz1_pre_trace = dfnz1.Signal_Pre_Timepoint_Averaged;
dfnz1_pre_sem = dfnz1.Signal_Pre_Timepoint_SEM;
dfnz1_post_trace = dfnz1.Signal_Post_Timepoint_Averaged;
dfnz1_post_sem = dfnz1.Signal_Post_Timepoint_SEM;
dfnz1_combined = [dfnz1_pre_trace dfnz1_post_trace];
dfnz1_sem = [dfnz1_pre_sem dfnz1_post_sem];
dfnz3 = combined_structure(6);
dfnz3_pre_trace = dfnz3.Signal_Pre_Timepoint_Averaged;
dfnz3_pre_sem = dfnz3.Signal_Pre_Timepoint_SEM;
dfnz3_post_trace = dfnz3.Signal_Post_Timepoint_Averaged;
dfnz3_post_sem = dfnz3.Signal_Post_Timepoint_SEM;
dfnz3_combined = [dfnz3_pre_trace dfnz3_post_trace];
dfnz3_sem = [dfnz3_pre_sem dfnz3_post_sem];
common_time = linspace(-10, 60, 105000); % create a common time vector for all traces with the same number of data points as the signal (negative time represents baseline so zero is the injection time)
figure;
hold on;
fig1 = plot(common_time, saline_combined, 'Color', [0 0 0], 'LineWidth', 2, 'DisplayName', 'Saline'); % plot each average trace
fill([common_time, fliplr(common_time)], [saline_combined - saline_sem, fliplr(saline_combined + saline_sem)], [0 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % shade SEM on either side of trace
fig2 = plot(common_time, fnz003_combined, 'Color', [1 0.627 0.251], 'LineWidth', 2, 'DisplayName', 'FNZ 0.03 mg/kg');
fill([common_time, fliplr(common_time)], [fnz003_combined - fnz003_sem, fliplr(fnz003_combined + fnz003_sem)], [1 0.627 0.251], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fig3 = plot(common_time, fnz01_combined, 'Color', [1 0.376 0], 'LineWidth', 2, 'DisplayName', 'FNZ 0.1 mg/kg');
fill([common_time, fliplr(common_time)], [fnz01_combined - fnz01_sem, fliplr(fnz01_combined + fnz01_sem)], [1 0.376 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fig4 = plot(common_time, dfnz03_combined, 'Color', [0.627 1 0.627], 'LineWidth', 2, 'DisplayName', 'DFNZ 0.3 mg/kg');
fill([common_time, fliplr(common_time)], [dfnz03_combined - dfnz03_sem, fliplr(dfnz03_combined + dfnz03_sem)], [0.627 1 0.627], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fig5 = plot(common_time, dfnz1_combined, 'Color', [0 1 0], 'LineWidth', 2, 'DisplayName', 'DFNZ 1 mg/kg');
fill([common_time, fliplr(common_time)], [dfnz1_combined - dfnz1_sem, fliplr(dfnz1_combined + dfnz1_sem)], [0 1 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fig6 = plot(common_time, dfnz3_combined, 'Color', [0 0.502 0], 'LineWidth', 2, 'DisplayName', 'DFNZ 3 mg/kg');
fill([common_time, fliplr(common_time)], [dfnz3_combined - dfnz3_sem, fliplr(dfnz3_combined + dfnz3_sem)], [0 0.502 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xline(0, '--'); % set a line at the injection time
plot([-10 -10], [-2 -3], 'LineWidth', 2, 'Color', 'k'); % plot a reference bar to show 1 z-score
plot([-10 -8], [-3 -3], 'LineWidth', 2, 'Color', 'k'); % plot a reference bar to show 2 minutes
axis off;
