function [data_norm] = normalizeFP3(input_data, siglog, type, plotlog, base_range, norm_range)
% This function uses a linear polynomial fit to normalize a sensor signal to a control signal to eliminate
% photobleaching and motion effects
% Input 1 (required): Input Data - A data structure containing fields Control_Signal,
% Sensor_Signal, Time_Vector, TTLs
% Input 2 (Required): Whole Signal Logical - Logical input tells the function whether
% to normalize based on the whole signal or a partial signal
% Input 3 (Required): Type - Character input that tells the function which
% type of output is desired (dff or zscore)
% Input 4 (Required): Plot Logical - Logical input to tell the function to
% plot the trace
% Input 5 (Optional): Baseline Range - Time range for baseline period
% Input 6 (Optional): Normalization Range - Time range of data to be
% normalized
% Output: data structure containing normalized trace, time vector, TTLs,
% and sampling rate

% Define inputs
x = inputParser;
addRequired(x, 'input_data', @isstruct);
addRequired(x, 'siglog', @islogical);
addRequired(x, 'type', @ischar);
addRequired(x, 'plotlog', @islogical);
addOptional(x, 'base_range', 0, @isvector);
addOptional(x, 'norm_range', 0, @isvector);

% Pull data out of input structure
sensor = input_data.Sensor_Signal;
control = input_data.Control_Signal;
time = input_data.Time_Vector;
ttls = input_data.TTLs;
sampling_rate = input_data.Sampling_Rate;

% Whole Signal Normalization
if siglog == 1
    % Use polyfit to define the parameters
    p = polyfit(control, sensor, 1);
    a = p(1); b = p(2);
    % Scale the control
    scaled_control = control * a + b;
    % Subtract scaled control and divde to correct for motion and
    % photobleaching, respectively, to get ΔF/F
    dff = (sensor - scaled_control) ./ scaled_control;
    % Calculate z-score
    med = median(dff);
    absolute_deviation = abs(dff - med);
    mad = median(absolute_deviation);
    z_score = (dff - median(dff)) / mad;
    % If specified, plot data
    % Plots z score with no baseline shift
    prompt1 = 'Do you want to shift the baseline (1 = yes, 0 = no)?';
    prompt2 = 'What value do you want to shift the baseline to?';
    shiftlog = input(prompt1);
    if shiftlog == 1
        shiftval = input(prompt2);
    end
    if strcmp(type, 'zscore') == 1 & shiftlog == 0
        timevec = time;
        signal = z_score;
    % Plots z score shifted so that the baseline is the value specified
    %by shiftval
    elseif strcmp(type, 'zscore') == 1 & shiftlog == 1
        baseline = z_score;
        mb = mean(baseline);
        diff1 = shiftval - mb;
        new_z = z_score + diff1;
        timevec = time;
        signal = new_z;
    % Plots ΔF/F shifted so that the baseline is the value
    % specified by shiftval
    elseif strcmp(type, 'dff') == 1 & shiftlog == 1
        baseline = dff;
        mb = mean(baseline);
        diff1 = shiftval - mb;
        new_dff = dff + diff1;
        timevec = time;
        signal = new_dff;
    % Plots the ΔF/F with no baseline shift
    elseif strcmp(type, 'dff') == 1 & shiftlog == 0
        timevec = time;
        signal = dff;
    end
end

% Normalization Based on a Baseline
if siglog == 0
% Ranges are assumed to be numeric vectors with the start and end time
% of the range
% The following section indexes the locations of the range starts and
% ends within the data
    base_start = base_range(1);
    base_start_log = time > base_start;
    post_base_start_locs = find(base_start_log == 1);
    base_start_loc = post_base_start_locs(1);
    base_end = base_range(2);
    base_end_log = time > base_end;
    post_base_end_locs = find(base_end_log == 1);
    base_end_loc = post_base_end_locs(1);
    norm_start = norm_range(1);
    norm_start_log = time >= norm_start;
    post_norm_start_locs = find(norm_start_log == 1);
    norm_start_loc = post_norm_start_locs(1);
    norm_end = norm_range(2);
    norm_end_log = time > norm_end;
    post_norm_end_locs = find(norm_end_log == 1);
    if ~ isempty(post_norm_end_locs)
        norm_end_loc = post_norm_end_locs(1);
    else
        norm_end_loc = length(time);
    end
    % Finds a polynomial fit only based on the baseline range specified
    p = polyfit(control(base_start_loc:base_end_loc), sensor(base_start_loc:base_end_loc), 1);
    a = p(1); b = p(2);
    % Applies the normalization only to the range specified
    scaled_control = control(norm_start_loc:norm_end_loc) * a + b;
    dff = (sensor(norm_start_loc:norm_end_loc) - scaled_control)./scaled_control;
    scaled_control_baseline = control(base_start_loc:base_end_loc) * a + b;
    baseline_dff = (sensor(base_start_loc:base_end_loc) - scaled_control_baseline)./scaled_control_baseline;
    % Calculates Z-Score for baseline
    medBASE = median(baseline_dff);
    absolute_deviation_base = abs(baseline_dff - medBASE);
    mad_base = median(absolute_deviation_base);
    z_score_base = (baseline_dff - median(baseline_dff)) / mad_base;
    % Calculates Z-Score for normalized signal
    med = median(dff);
    absolute_deviation = abs(dff - med);
    mad = median(absolute_deviation);
    z_score = (dff - median(dff)) / mad;
    % Ensures that the time vector matches the new normalized trace
    newtime = time(norm_start_loc:norm_end_loc);
    prompt1 = 'Do you want to shift the baseline (1 = yes, 0 = no)?';
    prompt2 = 'What value do you want to shift the baseline to?';
    shiftlog = input(prompt1);
    if shiftlog == 1
        shiftval = input(prompt2);
    end
    % Plots data if specified
    % Plots z score with no baseline shift
    if strcmp(type, 'zscore') == 1 & shiftlog == 0
        signal = z_score;
        timevec = newtime;
    % Plots z score shifted so that the baseline is the value specified
    %by shiftval
    elseif strcmp(type, 'zscore') == 1 & shiftlog == 1
        baseline = z_score(base_start_loc:base_end_loc);
        mb = mean(baseline);
        diff1 = shiftval - mb;
        new_z = z_score + diff1;
        signal = new_z;
        timevec = newtime;
    % Plots ΔF/F shifted so that the baseline is the value
    % specified by shiftval
    elseif strcmp(type,'dff') == 1 & shiftlog == 1
        baseline = dff(base_start_loc:base_end_loc);
        mb = mean(baseline);
        diff1 = shiftval - mb;
        new_dff = dff + diff1;
        signal = new_dff;
        timevec = newtime;
    % Plots the ΔF/F with no baseline shift
    elseif strcmp(type,'dff') == 1 & shiftlog == 0
        timevec = newtime;
        signal = dff;
    end
end
data_norm = struct('Signal_Vector', signal, 'Time_Vector', timevec, 'TTLs', {ttls}, 'Sampling_Rate', sampling_rate);
if plotlog == 1
    switch type
        case 'dff'
            figure;
            plot(timevec, signal);
            xlabel('Time (s)');
            ylabel('ΔF/F');
            title('ΔF/F of Normalized TDT Photometry Trace');
        case 'zscore'
            figure;
            plot(timevec, signal);
            xlabel('Time (s)');
            ylabel('Z-Score');
            title('Z-Score of Normalized TDT Photometry Trace');
    end
end
end
