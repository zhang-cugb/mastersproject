function [ structures] = importOPTVlog( boreholename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


filename = strcat(boreholename,'_structures.txt');
delimiter = '\t';
startRow = 3;

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Create output variable
structures = table(dataArray{1:end-1}, 'VariableNames', {'Depth','Azimuth','Dip','Aperture','Type'});

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

end

