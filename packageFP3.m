function [data_struct] = packageFP3(data_name, input_data, identifiers)
% Packagaes processed photometry data into a structure
% Input 1 (Required): Data Name - Name of the data structure to create/add to
% Input 2 (Required): Input Data - Data to be added to the structure
% Input 3 (Required): Identifiers - Logical input tells the function
% whether to ask for user input identifiers
% Output: Data structure containing relevant information in one row for
% each subject/event

% Define inputs as required
p = inputParser;
addRequired(p, 'data_name', @ischar);
addRequired(p, 'input_data', @isstruct);
addRequired(p, 'identifiers', @islogical);

fields = fieldnames(input_data);
tempcell = [];
if ~ exist(data_name, 'file')
    s = struct();
    structure_index = 1;
else
    tempStruct = load(data_name);
    varname = 's';
    fieldNames = fieldnames(tempStruct);
    eval([varname ' =tempStruct.' fieldNames{1} ';']);
    structure_size = width(s);
    structure_index = structure_size + 1;
end
if identifiers == 1
    prompt1 = 'What is the Animal ID?';
    prompt2 = 'What is the location?';
    prompt3 = 'What is the sensor?';
    prompt4 = 'What is the experiment?';
    s(structure_index).Animal_ID = input(prompt1);
    s(structure_index).Location = input(prompt2);
    s(structure_index).Sensor = input(prompt3);
    s(structure_index).Experiment = input(prompt4);
end
for a = 1:height(fields)
    field = fields{a};
    val = input_data.(field);
    s(structure_index).(field) = val;
end
data_struct = s;
