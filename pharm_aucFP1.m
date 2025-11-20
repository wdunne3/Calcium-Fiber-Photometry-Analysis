function auc_struct = pharm_aucFP1(input_data, window, max_post_window, start_ttl, plotlog)
% This function takes a normalized photometry trace with an associated ttl
% and creates a pre-ttl bucket and then divides the post ttl recordings
% into buckets, then takes the area under the curve

% Input 1: Input Data - structure containing Signal_Vector, Time_Vector,
% TTLs
% Input 2: Window - numeric input tells the function the winodw/bucket size
% in seconds
% Input 3: Max Post Window - numeric input tells the function how long
% after the timepoint to look at (in seconds)
% Input 4: Start TTL - character input tells the function which channel is
% the relevnt ttl
% Input 5: Plot Logical - logical input tells the function whether to plot
% data

% Output: Data structure containing the bucketed times, signals, and areas
% under the curve, as well as the bucket size used and relevant channel

% Define the inputs
p = inputParser;
addRequired(p, 'input_data', @isstruct);
addRequired(p, 'window', @isnumeric);
addRequired(p, 'max_post_window', @isnumeric);
addRequired(p, 'start_ttl', @ischar);
addOptional(p, 'plotlog', 0, @islogical);

% Pull signal information from the input data structure
signal = input_data.Signal_Vector;
time = input_data.Time_Vector;
ttls = input_data.TTLs;

% Use start ttl input to find the relevant timepoint
if ischar(start_ttl)
    num = str2num(start_ttl(3));
    idx = num + 1;
    pc = ttls{idx};
    timepoint = pc(1);
elseif isnumeric(start_ttl)
    timepoint = start_ttl;
end

% Finds the amount of time in seconds before the first timepoint
start = ceil(timepoint - time(1));

% Find the time in signal before and after timepoint, find auc after
% timepoint
max_post_time = max_post_window + timepoint;
post_tp_idx = find(time > timepoint & time <= max_post_time);
post_tp_times = time(post_tp_idx);
post_tp_signal = signal(post_tp_idx);
post_tp_auc = trapz(post_tp_times,post_tp_signal, 1);

% Create time bucket(s) for pre timepoint AUC analysis, find pre timepoint AUC
pre_num_bucs = floor(start / window);
pre_tp_idx = find(time < timepoint & time > timepoint - (window * pre_num_bucs));
pre_tp_times = time(pre_tp_idx);
pre_tp_signal = signal(pre_tp_idx);
pre_tp_auc = trapz(pre_tp_times,pre_tp_signal, 1);

% Create buckets for pre injection
min = window / 60;
aucbux = [];
sigbux = [];
timebux = [];
graphnames = [];
reverse_pre_num_bucs = [pre_num_bucs:-1:1];
for h = 1:pre_num_bucs
    pre_buc_idx = find(time <= timepoint - window * (h-1) & time > timepoint - window * (h));
    pre_bucket_time = time(pre_buc_idx); % add time, signal within each bucket
    pre_bucket = signal(pre_buc_idx);
    pre_auc_buc = trapz(pre_bucket_time, pre_bucket); % calculate the AUC for the bucket
    aucbux = [aucbux pre_auc_buc]; % add AUC, time, signal to vector for all buckets
    timebux = [timebux; pre_bucket_time];
    sigbux = [sigbux pre_bucket];
    name = [num2str(reverse_pre_num_bucs(h)*min) '-' num2str(reverse_pre_num_bucs(h)*min-min) ' min pre timepoint']; % create a bucket name
    graphname = {name};
    graphnames = [graphnames graphname]; % add bucket name to vector for all buckets
end

% Create custom post timepoint time buckets
num_bucs = max_post_window / window;
for i = 1:num_bucs
    buc_idx = find(time <= timepoint + window * i & time > timepoint + window * (i-1));
    bucket_time = time(buc_idx); % add time, signal within each bucket
    bucket = signal(buc_idx);
    auc_buc = trapz(bucket_time, bucket); % calculate AUC for each bucket
    aucbux = [aucbux auc_buc]; % add AUC, time, signal to vector for all buckets
    timebux = [timebux; bucket_time];
    sigbux = [sigbux bucket];
    name = [num2str(i*min - min) '-' num2str(i*min) ' min post timepoint']; % create a bucket name
    graphname = {name};
    graphnames = [graphnames graphname]; % add bucket name to vector for all buckets
end

% Add data to final output structure
auc_struct = struct('Bucket_Size', [num2str(min) ' minutes'], 'Post_Timepoint_Window', [num2str(max_post_window / 60) ' minutes'], 'Signal_Pre_Timepoint', pre_tp_signal, 'AUC_Pre_Timepoint', pre_tp_auc, 'Signal_Post_Timepoint', post_tp_signal, 'AUC_Post_Timepoint', post_tp_auc, 'Bucket_Labels', {graphnames}, 'Time_Buckets', timebux, 'Signal_Buckets', {sigbux}, 'AUC_Buckets', {aucbux}, 'Timepoint_Channel', start_ttl);
% Plot figures
if plotlog == 1
    % Plot signal with post timepoint AUC shaded
    figure;
    figure1 = plot(time, signal);
    hold on
    area(post_tp_times,post_tp_signal,0,"FaceColor",'r');
    figure1.Color = 'k';
    title('Normalized Photometry Trace');
    xlabel('Time (s)');
    ylabel('Signal');
    xline(timepoint, '--', 'Timepoint');
    % Plot bar graph of AUCs for each bucket
    figure;
    bar(graphnames,aucbux)
    title(['Area Under Curve For ' num2str(min) ' Minute Time Buckets']);
    xlabel('Bucket');
    ylabel('Area Under Curve');
end
