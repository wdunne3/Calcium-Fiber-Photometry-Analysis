%% Extract Raw Data From Structure, Downsample
idx = 1; % index data structure of raw data for individual entries
row = raw_structure(idx); % get one individual recording
time = row.Time_Vector; % pull data from individual entry
sensor = row.Sensor_Signal;
control = row.Control_Signal;
fs = row.Sampling_Rate; 

% Define parameters for baseline length, recording length, window length
baseline_window = 600; 
post_injection_window = 3600;
bucket_window = 300;

% Downsample to 20 Hz
target_fs = 20;
if fs > target_fs
    downFactor = floor(fs/target_fs);
    sensor_ds = downsample(sensor, downFactor);
    control_ds = downsample(control, downFactor);
    time = downsample(time, downFactor);
    fs = fs/downFactor;
end

% Extract TTL (PC1)
ttl_mat = row.TTLs;
pc1 = ttl_mat{2};
drug_event = pc1(1);
time = time - drug_event; % adjust so that injection is 0
% Baseline: Specified time before drug injection
baseline_period = time < 0 & time >= -baseline_window;

% ΔF/F calculation
x_mean = mean(control_ds(baseline_period));
x_std = std(control_ds(baseline_period));
scaled_isosbestic = (control_ds(baseline_period)-x_mean)/x_std;

try
    p = polyfit(scaled_isosbestic, sensor_ds(baseline_period), 1);
catch
    warning('Polynomial fit failed, using default scaling.');
    p = [0 0];
end


fitted_control = polyval(p, (control_ds-x_mean)/x_std); % apply polynomial fit to control
dF_F = (sensor_ds - fitted_control)./fitted_control; % subtract and divide fitted control to get dff
figure;
plot(time, dF_F); % visualize dff

% Iterative photobleaching correction (for spike detection)
detrended_signal = zeros(size(dF_F));
window_size = 15;
num_iterations = 100;
for iter = 1:num_iterations
    offset = mod(iter-1,window_size)+1;
    mov_avg = movmean(dF_F,window_size,'Endpoints','shrink');
    normalized_window = dF_F - [NaN(1,offset-1), mov_avg(1:end-offset+1)];
    detrended_signal = detrended_signal + normalized_window/num_iterations;
end


figure;
plot(time, detrended_signal); % visualize detrended signal
xline(0, '--'); % plot lines to show injection, beginng, end of period of interest
xline(-baseline_window);
xline(post_injection_window);
% Post-injection
post_drug_period = time >= 0 & time < post_injection_window; % mask to find one hour after injection time
drug_dF_F = detrended_signal(post_drug_period);
drug_time = time(post_drug_period);

% Pre-injecton
pre_drug_period = time <= 0 & time >  -baseline_window; % mask to find baseline period before injection
baseline_dF_F = detrended_signal(pre_drug_period);
baseline_time = time(pre_drug_period);

% Spike detection (high-threshold)
hold_time = 1;
drug_thresh = mean(drug_dF_F)+3*std(drug_dF_F); % set the detection threshold for spikes > 3 std above the mean signal
drug_spike_idx = find(drug_dF_F > drug_thresh); % index spikes identified, find associated times
drug_spike_times = drug_time(drug_spike_idx);
baseline_thresh = mean(baseline_dF_F) + 3 * std(baseline_dF_F); % detect spikes in the baseline
baseline_spike_idx = find(baseline_dF_F > baseline_thresh);
baseline_spike_times = baseline_time(baseline_spike_idx);

if ~isempty(drug_spike_times) % filter spikes
    drug_spike_times_filtered = drug_spike_times(1);
    for i=2:length(drug_spike_times)
        if drug_spike_times(i)-drug_spike_times_filtered(end) >= hold_time
            drug_spike_times_filtered = [drug_spike_times_filtered; drug_spike_times(i)];
        end
    end
else
    drug_spike_times_filtered = [];
end

if ~isempty(baseline_spike_times) % filter baseline spikes
    baseline_spike_times_filtered = baseline_spike_times(1);
    for i=2:length(baseline_spike_times)
        if baseline_spike_times(i)-baseline_spike_times_filtered(end) >= hold_time
            baseline_spike_times_filtered = [baseline_spike_times_filtered; baseline_spike_times(i)];
        end
    end
else
    baseline_spike_times_filtered = [];
end

% Spike amplitude and width
drug_amp = zeros(length(drug_spike_times_filtered),1);
drug_width = zeros(length(drug_spike_times_filtered),1);
for i=1:length(drug_spike_times_filtered)
    t = drug_spike_times_filtered(i);
    win = (drug_time >= t-1 & drug_time <= t+1); % set a one-second window around each spike
    [drug_amp(i),~] = max(drug_dF_F(win));
    w = find(drug_dF_F(win) > drug_thresh, 1, 'last') - find(drug_dF_F(win) > drug_thresh, 1, 'first');
    drug_width(i) = w/fs;
end

% Amplitude and Width for Baseline
baseline_amp = zeros(length(baseline_spike_times_filtered),1);
baseline_width = zeros(length(baseline_spike_times_filtered),1);
for i=1:length(baseline_spike_times_filtered)
    t = baseline_spike_times_filtered(i);
    win = (baseline_time >= t-1 & baseline_time <= t+1);
    [baseline_amp(i),~] = max(baseline_dF_F(win));
    w = find(baseline_dF_F(win) > baseline_thresh, 1, 'last') - find(baseline_dF_F(win) > baseline_thresh, 1, 'first');
    baseline_width(i) = w/fs;
end

baseline_rate = numel(baseline_spike_times_filtered) / baseline_window; % divide by baseline length to find spike rate in Hz
drug_rate = numel(drug_spike_times_filtered) / post_injection_window; % divide by treatment length to find spike rate in Hz


% ΔF/F Event analysis (lower-threshold transients)
event_thresh = mean(drug_dF_F) + 1.5*std(drug_dF_F); % sets lower threshold at 1.5 std above mean
above_thresh = drug_dF_F > event_thresh;
above_thresh = above_thresh(:);  % Ensure column vector

diff_var = [0; above_thresh; 0]; 
d_event = diff(diff_var);
event_starts = find(d_event==1); % find start and end indices and times for events
event_ends = find(d_event==-1)-1;
num_events = length(event_starts);
event_start_times = drug_time(event_starts);
event_start_times = event_start_times(:); % ensure column
event_end_times = drug_time(event_ends);
event_end_times = event_end_times(:); % ensure column
event_windows = [event_start_times event_end_times];

event_amp = zeros(num_events,1);
event_width = zeros(num_events,1);
event_area = zeros(num_events,1);
event_times = zeros(num_events,1); % add exact times for max amp

dt = 1/fs; % sampling interval
for i=1:num_events
    idx = event_starts(i):event_ends(i);
    event_amp(i) = max(drug_dF_F(idx));
    event_idc = find(drug_dF_F == event_amp(i)); % find exact index for max amp
    event_times(i) = drug_time(event_idc); % add exact time for max amp
    event_width(i) = length(idx)*dt;
    event_area(i) = sum(drug_dF_F(idx))*dt;  % safe area calculation
end
event_rate = num_events / post_injection_window; % get frequency in Hz (divide by drug time)

% Repeat event analysis for baseline
baseline_event_thresh = mean(baseline_dF_F) + 1.5 * std(baseline_dF_F);
baseline_above_thresh = baseline_dF_F > baseline_event_thresh;
baseline_above_thresh = baseline_above_thresh(:);

baseline_diff_var = [0; baseline_above_thresh; 0];
baseline_d_event = diff(baseline_diff_var);
baseline_event_starts = find(baseline_d_event==1);
baseline_event_ends = find(baseline_d_event==-1)-1;
baseline_num_events = length(baseline_event_starts);
baseline_event_start_times = baseline_time(baseline_event_starts);
baseline_event_start_times = baseline_event_start_times(:);
baseline_event_end_times = baseline_time(baseline_event_ends);
baseline_event_end_times = baseline_event_end_times(:);
baseline_event_windows = [baseline_event_start_times baseline_event_end_times];

baseline_event_amp = zeros(baseline_num_events,1);
baseline_event_width = zeros(baseline_num_events,1);
baseline_event_area = zeros(baseline_num_events,1);
baseline_event_times = zeros(baseline_num_events, 1); % add exact times for max amp

baseline_dt = 1/fs; % sampling interval
for i=1:baseline_num_events
    idx = baseline_event_starts(i):baseline_event_ends(i);
    baseline_event_amp(i) = max(baseline_dF_F(idx));
    baseline_event_idc = find(baseline_dF_F == baseline_event_amp(i)); % find exact index for max amp
    baseline_event_times(i) = baseline_time(baseline_event_idc); % add exact time for max amp
    baseline_event_width(i) = length(idx)*baseline_dt;
    baseline_event_area(i) = sum(baseline_dF_F(idx))*baseline_dt;  % safe area calculation
end
baseline_event_rate = baseline_num_events / baseline_window; % get frequency in Hz (divide by baseline time)

% Separate Spike and Events into Specified Blocks
graphnames = []; % create empty vectors to place bucketed data
bucket_signals = [];
bucket_amps = [];
av_bucket_amps = [];
bucket_widths = [];
av_bucket_widths = [];
bucket_rates = [];
bucket_event_amps = [];
av_bucket_event_amps = [];
bucket_event_widths = [];
av_bucket_event_widths = [];
bucket_event_areas = [];
av_bucket_event_areas = [];
bucket_event_rates = [];
% Make buckets for baseline
baseline_num_buckets = baseline_window / bucket_window; % number of buckets in baseline period
for a = [baseline_num_buckets:-1:1]
    low = -bucket_window * a; % set range of times for buckets
    high = -bucket_window * (a-1);
    baseline_bucket_idx = find(time >= low & time < high); % mask to find indices
    baseline_bucket_time = time(baseline_bucket_idx);
    baseline_bucket_signal = detrended_signal(baseline_bucket_idx);
    bucket_signals = [bucket_signals; baseline_bucket_signal]; % add signal to bucket
    baseline_spike_mask = baseline_spike_times_filtered >= baseline_bucket_time(1) & baseline_spike_times_filtered <= baseline_bucket_time(end); % mask to find spikes within time bucket
    baseline_bucket_amp = baseline_amp(baseline_spike_mask);
    av_bucket_amp = mean(baseline_bucket_amp); % find average spike amplitude for time bucket
    bucket_amps = [bucket_amps; {baseline_bucket_amp}]; % add all bucket spikes
    av_bucket_amps = [av_bucket_amps; av_bucket_amp]; % add average value for bucket
    baseline_bucket_width = baseline_width(baseline_spike_mask);
    av_bucket_width = mean(baseline_bucket_width); % repeat for spike measures
    bucket_widths = [bucket_widths; {baseline_bucket_width}];
    av_bucket_widths = [av_bucket_widths; av_bucket_width];
    baseline_bucket_rate = numel(find(baseline_spike_mask == 1)) / bucket_window; % calculate spike rate within bucket
    bucket_rates = [bucket_rates; baseline_bucket_rate];
    baseline_event_mask = baseline_event_times >= baseline_bucket_time(1) & baseline_event_times <= baseline_bucket_time(end); % mask to find events within bucket
    bucket_event_amp = baseline_event_amp(baseline_event_mask);
    av_bucket_event_amp = mean(bucket_event_amp); % repeat average values for events
    bucket_event_amps = [bucket_event_amps; {bucket_event_amp}];
    av_bucket_event_amps = [av_bucket_event_amps; av_bucket_event_amp];
    bucket_event_width = baseline_event_width(baseline_event_mask);
    av_bucket_event_width = mean(bucket_event_width);
    bucket_event_widths = [bucket_event_widths; {bucket_event_width}];
    av_bucket_event_widths = [av_bucket_event_widths; av_bucket_event_width];
    bucket_event_area = baseline_event_area(baseline_event_mask);
    av_bucket_event_area = mean(bucket_event_area);
    bucket_event_areas = [bucket_event_areas; {bucket_event_area}];
    av_bucket_event_areas = [av_bucket_event_areas; av_bucket_event_area];
    bucket_event_rate = numel(bucket_event_amp) / bucket_window; % find rate of events within bucket
    bucket_event_rates = [bucket_event_rates; bucket_event_rate];
    name = [num2str(bucket_window*a/60) '-' num2str(bucket_window*(a-1)/60) ' min pre timepoint']; % create neame to plot on a bar graph
    graphname = {name};
    graphnames = [graphnames graphname];
end

drug_num_buckets = post_injection_window / bucket_window; % number of buckets within post-injection window
for c = 1:drug_num_buckets
    low = bucket_window * (c-1); % set range of times within bucket
    high = bucket_window * c;
    drug_bucket_idx = find(time >= low & time < high); % mask to find bucket
    drug_bucket_time = time(drug_bucket_idx);
    drug_bucket_signal = detrended_signal(drug_bucket_idx);
    bucket_signals = [bucket_signals; drug_bucket_signal];
    drug_spike_mask = drug_spike_times_filtered >= drug_bucket_time(1) & drug_spike_times_filtered <= drug_bucket_time(end); % mask to find spikes within the bucket
    if isempty(find(drug_spike_mask == 1, 1)) % in case there are no spikes within a bucket
        drug_bucket_amp = 0;
        av_bucket_amp = 0;
        drug_bucket_width = 0;
        av_bucket_width = 0;
        drug_bucket_rate = 0;
        bucket_amps = [bucket_amps; drug_bucket_amp];
        av_bucket_amps = [av_bucket_amps; av_bucket_amp];
        bucket_widths = [bucket_widths; drug_bucket_width];
        av_bucket_widths = [av_bucket_widths; av_bucket_width];
        bucket_rates = [bucket_rates; drug_bucket_rate];
    else
        drug_bucket_amp = drug_amp(drug_spike_mask);
        av_bucket_amp = mean(drug_bucket_amp); % find spike measures, averages within bucket
        bucket_amps = [bucket_amps; drug_bucket_amp];
        av_bucket_amps = [av_bucket_amps; av_bucket_amp];
        drug_bucket_width = drug_width(drug_spike_mask);
        av_bucket_width = mean(drug_bucket_width);
        bucket_widths = [bucket_widths; drug_bucket_width];
        av_bucket_widths = [av_bucket_widths; av_bucket_width];
        drug_bucket_rate = numel(find(drug_spike_mask == 1)) / bucket_window; % calculate spike rate within bucket
        bucket_rates = [bucket_rates; drug_bucket_rate];
    end
    bucket_event_mask = event_times >= drug_bucket_time(1) & event_times <= drug_bucket_time(end); % mask to find events within bucket
    if isempty(find(bucket_event_mask == 1, 1)) % in case there are not events within the window
        bucket_event_amp = 0;
        av_bucket_event_amp = 0;
        bucket_event_amps = [bucket_event_amps; {bucket_event_amp}];
        av_bucket_event_amps = [av_bucket_event_amps; av_bucket_event_amp];
        bucket_event_width = 0;
        av_bucket_event_width = 0;
        bucket_event_widths = [bucket_event_widths; {bucket_event_width}];
        av_bucket_event_widths = [av_bucket_event_widths; av_bucket_event_width];
        bucket_event_area = 0;
        av_bucket_event_area = 0;
        bucket_event_areas = [bucket_event_areas; {bucket_event_area}];
        av_bucket_event_areas = [av_bucket_event_areas; av_bucket_event_area];
        bucket_event_rate = 0;
        bucket_event_rates = [bucket_event_rates; bucket_event_rate];
    else
        bucket_event_amp = event_amp(bucket_event_mask);
        av_bucket_event_amp = mean(bucket_event_amp); % find event measures and averages
        bucket_event_amps = [bucket_event_amps; {bucket_event_amp}];
        av_bucket_event_amps = [av_bucket_event_amps; av_bucket_event_amp];
        bucket_event_width = event_width(bucket_event_mask);
        av_bucket_event_width = mean(bucket_event_width);
        bucket_event_widths = [bucket_event_widths; {bucket_event_width}];
        av_bucket_event_widths = [av_bucket_event_widths; av_bucket_event_width];
        bucket_event_area = event_area(bucket_event_mask);
        av_bucket_event_area = mean(bucket_event_area);
        bucket_event_areas = [bucket_event_areas; {bucket_event_area}];
        av_bucket_event_areas = [av_bucket_event_areas; av_bucket_event_area];
        bucket_event_rate = numel(bucket_event_amp) / bucket_window; % find event rate within bucket
        bucket_event_rates = [bucket_event_rates; bucket_event_rate];
    end
    name = [num2str(bucket_window*(c-1)) '-' num2str(300*c) ' min post timepoint']; % create name for plotting bucket
    graphname = {name};
    graphnames = [graphnames graphname];
end

% Normalize values as a percent of baseline
average_amp_baseline = mean(baseline_amp); % find averages of baseline spike and event measure to normalize buckets
average_width_baseline = mean(baseline_width);
av_bucket_widths_perc = av_bucket_widths / average_width_baseline * 100; % adjust to percent baseline
av_bucket_amps_perc = av_bucket_amps / average_amp_baseline * 100; % adjust to percent baseline
bucket_rates_perc = bucket_rates / baseline_rate * 100; % adjust to percent baseline
bucket_event_rates_perc = bucket_event_rates / baseline_event_rate * 100; % adjust to percent baseline
average_baseline_event_amp = mean(baseline_event_amp);
average_baseline_event_width = mean(baseline_event_width);
average_baseline_event_area = mean(baseline_event_area);
av_bucket_event_amps_perc = av_bucket_event_amps / average_baseline_event_amp * 100; % adjust to percent baseline
av_bucket_event_widths_perc = av_bucket_event_widths / average_baseline_event_width * 100; % adjust to percent baseline
av_bucket_event_areas_perc = av_bucket_event_areas / average_baseline_event_area * 100; % adjust to percent baseline
average_baseline_event_amp = mean(baseline_event_amp); % average amplitude, width, area of events for normalization to baseline
average_baseline_event_width = mean(baseline_event_width);
average_baseline_event_area = mean(baseline_event_area);
baseline_event_amp_perc = baseline_event_amp / average_baseline_event_amp * 100; % convert to percent baseline
event_amp_perc = event_amp / average_baseline_event_amp * 100; % convert to percent baseline
baseline_event_width_perc = baseline_event_area / average_baseline_event_area * 100; % convert to percent baseline
event_width_perc = event_area / average_baseline_event_area * 100; % convert to percent baseline
baseline_event_area_perc = baseline_event_area / average_baseline_event_area * 100; % convert to percent baseline
event_area_perc = event_area / average_baseline_event_area * 100; % convert to percent baseline

% Create a structure outputting all values
spike_event_info = struct('Baseline_Signal', baseline_dF_F, 'Treatment_Signal', drug_dF_F, 'Bucket_Labels', {graphnames}); % create a data structure with spike and event info
spike_event_info.Baseline_Spike_Rate_Hz = baseline_rate; % add spike and event info to structure
spike_event_info.Treatment_Spike_Rate_Hz = drug_rate;
spike_event_info.Baseline_Event_Rate_Hz = baseline_event_rate;
spike_event_info.Treatment_Event_Rate_Hz = event_rate;
spike_event_info.Bucket_Signals = bucket_signals;
spike_event_info.Bucket_Spike_Amps = bucket_amps;
spike_event_info.Average_Bucket_Spike_Amps = av_bucket_amps;
spike_event_info.Average_Bucket_Spike_Amps_Percent = av_bucket_amps_perc;
spike_event_info.Bucket_Spike_Widths = bucket_widths;
spike_event_info.Average_Bucket_Spike_Widths = av_bucket_widths;
spike_event_info.Average_Bucket_Spike_Widths_Percent = av_bucket_widths_perc;
spike_event_info.Bucket_Spike_Rates = bucket_rates;
spike_event_info.Bucket_Spike_Rates_Percent = bucket_rates_perc;
spike_event_info.Bucket_Event_Amps = bucket_event_amps;
spike_event_info.Average_Bucket_Event_Amps = av_bucket_event_amps;
spike_event_info.Average_Bucket_Event_Amps_Percent = av_bucket_event_amps_perc;
spike_event_info.Bucket_Event_Widths = bucket_event_widths;
spike_event_info.Average_Bucket_Event_Widths = av_bucket_event_widths;
spike_event_info.Average_Bucket_Event_Widths_Percent = av_bucket_event_widths_perc;
spike_event_info.Bucket_Event_Areas = bucket_event_areas;
spike_event_info.Average_Bucket_Event_Areas = av_bucket_event_areas;
spike_event_info.Average_Bucket_Event_Areas_Percent = av_bucket_event_areas_perc;
spike_event_info.Bucket_Event_Rates = bucket_event_rates;
spike_event_info.Bucket_Event_Rates_Percent = bucket_event_rates_perc;
%% Visualize Data (Using Spike Times and Amps)
full_time = [baseline_time drug_time]; % concatenate baseline and post-drug times
full_dF_F = [baseline_dF_F drug_dF_F]; % concatenate baseline and post-drug signals
figure;
subplot(2,1,1); % create two plots to visualize spikes and events
plot(full_time, full_dF_F);
hold on;
if ~isempty(baseline_spike_times_filtered)
    scatter(baseline_spike_times_filtered, baseline_amp, 'g','filled'); % plot baseline spikes as markers
end
if ~isempty(drug_spike_times_filtered)
    scatter(drug_spike_times_filtered, drug_amp, 'r','filled'); % plot post-drug spikes as markers
end
xlabel('Time (min)'); 
xline(0, '--', 'Injection'); % set a line to show where injection happened
ylabel('ΔF/F (detrended)');
title('Detrended ΔF/F with Spikes (Baseline: Green, Treatment: Red)');
grid on;
subplot(2,1,2); % create a second plot for events
hold on;
plot(full_time, full_dF_F);
% Plot events (lower-threshold)
if ~isempty(event_starts)
    % Map event_starts to indices in full_dF_F_raw
    scatter(baseline_event_times, baseline_event_amp, 'o','MarkerEdgeColor',[1 0.5 0],'MarkerFaceColor',[1 0.5 0]); % plot baseline events as markers
end
if ~isempty(baseline_event_starts)
    % Map event_starts to indices in full_dF_F_raw
    scatter(event_times, event_amp, 'o','MarkerEdgeColor',[1 1 0],'MarkerFaceColor',[1 1 0]); % plot post-drug events as markers
end
xlabel('Time (min)');
xline(0, '--', 'Injection'); % Set a line to show where injection happened
ylabel('ΔF/F (detrended)');
title('Detrended ΔF/F with Events (Baseline: Orange, Treatment: Yellow)');
grid on;

%% Save Data and Add to Structure
close all;
spikes_events = packageFP3('gcampopioidspikes.mat', spike_event_info, 1); % save individual recording spike data to larger structure
prompt = 'File Save Name: ';
save_name = input(prompt); % sanity check to ensure the correct file name is used
save(save_name, "spikes_events");
