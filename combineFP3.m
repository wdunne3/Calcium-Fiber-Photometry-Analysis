function [averaged] = combineFP3(input_data, identifiers)
% This function receives a full data structure of representative data from
% an experiment and averages it for all subjects with the same sensor in
% the same location

% Input 1 (Required): Input Data - data structure of trials from an
% experiment. Can be added as either a structure or a file name
% corresponding to a saved structure
% Input 2 (Required): Identifiers - nx2 cell matrix in which the left
% column contains a character vector for the field of interest and the
% right column contains the value of interest (ex: [{'Sensor'} {'DA'};
% {'Location'} {'DS'}]

% Output: Data structure containing non-numerical fields from input data
% and averages and SEM for numerical fields

% Define inputs
p = inputParser;
addRequired(p, 'input_data', @isstruct);
addRequired(p, 'analysis', @ischar);
addRequired(p, 'identifiers', @ischar);
addOptional(p, 'window', 0, @isnumeric);

% If entered as a filename, loads the input data structure
if ischar(input_data) == 1
   loaded = load(input_data);
   fieldcell = fieldnames(loaded);
   fieldname = char(fieldcell);
   input_data = loaded.(fieldname);
end

% Finds all fields in input, creates blank vectors to sort
fields = fieldnames(input_data);
numfields = length(fields);
blankmat = [];
idxvec = [];

% Iterates through input structure to find all rows with matching
% identifiers
n = height(identifiers);
for z = 1:width(input_data)
    row = input_data(z);
    compvec = [];
    for m = 1:n
        id = identifiers(m, :);
        field = id{1};
        val = id{2};
        log = strcmp(row.(field), val);
        compvec = [compvec; log];
    end
    if isempty(find(compvec == false))
        last_idx = z;
        field_vector = [];
        idxs = [];
        char_idx = [];
        for c = 1:numfields
            fieldname = fields{c};
            fieldval = row.(fieldname);
            blankvec = [];
            if ~ischar(fieldval) 
                if iscolumn(fieldval)
                    fieldval = fieldval';
                end
                blankvec = {fieldval};
                field_vector = [field_vector blankvec];
                idxs = [idxs c];
            else
                char_idx = [char_idx c];
            end
        end
        blankmat = [blankmat; field_vector];
    end
end
names = fields';
tablenames = names(idxs);
charnames = names(char_idx);
table = cell2table(blankmat);
table.Properties.VariableNames = tablenames;
last = input_data(last_idx);
for b = 1:length(charnames)
    na = charnames{b};
    averaged.(na) = last.(na);
end
for g = 1:n
    tv = identifiers(g, :);
    field = tv{1};
    val = tv{2};
    averaged.(field) = val;
end

% Separates out numerical and non-numerical fields, finds mean and SEM for
% numerical fields
for h = 1:width(table)
    data = table{:, h};
    num_subjects = height(table);
    if height(data) == 1
        data_mean = data;
        data_sem = std(data, 0, 1) / sqrt(num_subjects);
        header = table.Properties.VariableNames{h};
        header2 = [header '_Averaged'];
        header3 = [header '_SEM'];
        averaged.(header) = data;
        averaged.(header2) = data_mean;
        averaged.(header3) = data_sem;
    elseif ~isnumeric(data)
        header = table.Properties.VariableNames{h};
        averaged.(header) = data;
    else
        data_mean = mean(data);
        data_sem = std(data, 0, 1) / sqrt(num_subjects);
        header = table.Properties.VariableNames{h};
        header2 = [header '_Averaged'];
        header3 = [header '_SEM'];
        averaged.(header) = data;
        averaged.(header2) = data_mean;
        averaged.(header3) = data_sem;
    end 
end
end
