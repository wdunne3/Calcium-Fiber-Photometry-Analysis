function [rawdata] = rawTDT1(filepath, beginning_trim, end_trim, plotlog, ttls_to_plot)
% Plot Raw TDT Photometry Data with TTLs
% Input 1 (required): Filepath - folder where the raw data is stored
% Input 2 (required): Beginning Trim - time to trim off of the beginning of
% the recording
% Input 3 (required): End Trim - time to trim off the end of the recording
% Input 4 (optional): Plot Logical - a 0 or 1 that indicates whether or not
% to plot raw data
% Input 5 (optional): TTLs to plot - if plotting data, a character vector that lists each
% desired TTL channel (PC0, PC1, etc) separated by commas
% Output: data structure containing control signal, sensor signal, time
% vector, sampling rate, and TTLs

% Define inputs as required and optional
p = inputParser;
addRequired(p, 'filepath', @ischar);
addRequired(p, 'beginning_trim', @isnumeric);
addRequired(p, 'end_trim', @isnumeric);
addOptional(p, 'plotlog', 0, @islogical);
addOptional(p, 'ttls_to_plot', [], @ischar);

% Remove double quotes from filename, Read TDT data into a structure
quotes = find(filepath == '"');
filepath(quotes) = [];
filepath = char(filepath);
data = TDTbin2mat(filepath);

% Extract data streams from A or C channels, create a time vector to
% match the raw photometry data
if isfield(data.streams, 'x465C')
    sensor_raw = double(data.streams.x465C.data);
    control_raw = double(data.streams.x405C.data);
    sampling_rate = round(data.streams.x405C.fs);
    times = 1:length(sensor_raw);
    times = times/sampling_rate+data.streams.x405C.startTime;
elseif isfield(data.streams, 'x465A')
    sensor_raw = double(data.streams.x465A.data);
    control_raw = double(data.streams.x405A.data);
    sampling_rate = round(data.streams.x405A.fs);
    times = 1:length(sensor_raw);
    times = times/sampling_rate+data.streams.x405A.startTime;
end

% Add existing TTLs to a matrix
ttl_matrix = [];
if isfield(data.epocs, 'PC0_')
    pc0_ttls = data.epocs.PC0_.onset;
    ttl_matrix = [ttl_matrix; {pc0_ttls}];
else 
    ttl_matrix = [ttl_matrix; {'No PC0 TTLs'}];
end
if isfield(data.epocs, 'PC1_')
    pc1_ttls = data.epocs.PC1_.onset;
    ttl_matrix = [ttl_matrix; {pc1_ttls}];
else
    ttl_matrix = [ttl_matrix; {'No PC1 TTLs'}];
end
if isfield(data.epocs, 'PC2_')
    pc2_ttls = data.epocs.PC2_.onset;
    ttl_matrix = [ttl_matrix; {pc2_ttls}];
else
    ttl_matrix = [ttl_matrix; {'No PC2 TTLs'}];
end
if isfield(data.epocs, 'PC3_')
    pc3_ttls = data.epocs.PC3_.onset;
    ttl_matrix = [ttl_matrix; {pc3_ttls}];
else
    ttl_matrix = [ttl_matrix; {'No PC3 TTLs'}];
end

% Trim beginning and end of recording
bt = beginning_trim * sampling_rate;
et = end_trim * sampling_rate;
sensor_signal = sensor_raw(bt+1:length(times)-et);
control_signal = control_raw(bt+1:length(times)-et);
t = times(bt+1:length(times)-et);

% Generate a structure containing all relevant raw data
rawdata = struct('Control_Signal', control_signal, 'Sensor_Signal', sensor_signal, 'Time_Vector', t, 'Sampling_Rate', sampling_rate, 'TTLs', {ttl_matrix}, 'Raw_Data_File', filepath);

% Plot raw data with TTLs specified (if plotting)
if exist("plotlog", 'var')
    if plotlog == 1
        figure;
        hold on;
        line_sensor = plot(t, sensor_signal, 'g', 'DisplayName', 'Sensor');
        line_control = plot(t, control_signal, 'm', 'DisplayName', 'Control');
        labels = [line_sensor line_control];

        % Find which TTLs are specified, plot them as lines
        if ~ isempty(strfind(ttls_to_plot, 'PC0')) & exist("pc0_ttls", 'var')
            x0 = xline(pc0_ttls(1), '--', 'DisplayName', 'PC0');
            labels = [labels x0];

            % Ensures that the display name is only in the legend once, not
            % for every instance of the TTL
            if length(pc0_ttls) > 1
                xline(pc0_ttls(2:end), '--', 'HandleVisibility', 'off');
            end
        end
        if ~ isempty(strfind(ttls_to_plot, 'PC1')) & exist("pc1_ttls", 'var')
            x1 = xline(pc1_ttls(1), '-', 'DisplayName', 'PC1');
            labels = [labels x1];
            if length(pc1_ttls) > 1
                xline(pc1_ttls(2:end), '-', 'HandleVisibility', 'off');
            end
        end
        if ~ isempty(strfind(ttls_to_plot, 'PC2')) & exist("pc2_ttls", 'var')
            x2 = xline(pc2_ttls(1), ':', 'DisplayName', 'PC2');
            labels = [labels x2];
            if length(pc2_ttls) > 1
                xline(pc2_ttls(2:end), ':', 'HandleVisibility', 'off')
            end
        end
        if ~ isempty(strfind(ttls_to_plot, 'PC3')) & exist("pc3_ttls", 'var')
            x3 = xline(pc3_ttls(1), '-.', 'DisplayName', 'PC3');
            labels = [labels x3];
            if length(pc3_ttls) > 1
                xline(pc3_ttls(2:end), '-.', 'HandleVisibility', 'off')
            end
        end
        title('Raw Data');
        xlabel('Seconds');
        ylabel('Signal (mV)');
        legend(labels);
    end
end
end
