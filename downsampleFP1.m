function [data_ds] = downsampleFP1(input_data, target_Hz, plotlog)
% Downsamples a dataset to the specified target frequency
% Input 1 (Required): Structure of data to be downsampled. Must contain
% Control_Signal, Sensor_Signal, Time_Vector, and Sampling_Rate fields
% Input 2 (Required): Target frequency for final data, put in as a number
% Input 3 (Optional): Plot logical. Logical input determines whether the function
% generates a plot showing the control and sensor traces
% Output: Data structure containing downsampled signal, control, and time
% vectors, sampling rate, TTLs, and raw file

% Define required inputs
p = inputParser;
addRequired(p, 'input_data', @isstruct);
addRequired(p, 'target_Hz', @isnumeric);
addOptional(p, 'plotlog', false, @islogical);

% Pull data from input structure
control = input_data.Control_Signal;
sensor = input_data.Sensor_Signal;
time = input_data.Time_Vector;
sampling_rate = input_data.Sampling_Rate;
ttls = input_data.TTLs;
file = input_data.Raw_Data_File;

% Extend dataset to eliminate end effects of resample function
if iscolumn(time)
    time = time';
end
time_ext = [(-20:-1)*target_Hz time (1:20)*target_Hz+time(end)];
signal_ext = (interp1(time,sensor.',time_ext,'linear','extrap')).';
control_ext = (interp1(time,control.',time_ext,'linear','extrap')).';

% Use resample function to downsample to desired frequency
aligned_sensor = resample(signal_ext, target_Hz, sampling_rate);
aligned_control = resample(control_ext, target_Hz, sampling_rate);
aligned_time = resample(time_ext, target_Hz, sampling_rate);

% Trim extended data back to original dataset
aligned_sensor(1:20) = [];
aligned_sensor(end-20:end) = [];
aligned_control(1:20) = [];
aligned_control(end-20:end) = [];
aligned_time(1:20) = [];
aligned_time(end-20:end) = [];

% Package downsampled data into a structure for output
data_ds = struct('Control_Signal', aligned_control, 'Sensor_Signal', aligned_sensor, 'Time_Vector', aligned_time, 'Sampling_Rate', target_Hz, 'TTLs', {ttls}, 'Raw_Data_File', file);

% Plot
if ~ exist("plotlog", 'var')
    plotlog = false;
end
if plotlog == 1
    figure;
    hold on;
    plot(aligned_time, aligned_sensor, 'g', 'LineWidth', 2);
    plot(aligned_time, aligned_control, 'm', 'LineWidth', 2);
end
end
