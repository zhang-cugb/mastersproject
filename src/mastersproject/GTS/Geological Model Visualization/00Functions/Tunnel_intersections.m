function [ SZ_tunnel] = Tunnel_intersections(GTS_coordinates)
filename = 'Tunnel_intersections.txt';
delimiter = '\t';
formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,7]
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end
rawNumericColumns = raw(:, [1,2,3,4,5,7]);
rawStringColumns = string(raw(:, 6));
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";

% Create output variable
SZ_tunnel = table;
SZ_tunnel.x = cell2mat(rawNumericColumns(:, 1));
SZ_tunnel.y= cell2mat(rawNumericColumns(:, 2));
SZ_tunnel.z = cell2mat(rawNumericColumns(:, 3));
SZ_tunnel.Azimuth = cell2mat(rawNumericColumns(:, 4));
SZ_tunnel.Dip = cell2mat(rawNumericColumns(:, 5));
SZ_tunnel.Tunnel = categorical(rawStringColumns(:, 1));
SZ_tunnel.SZ = cell2mat(rawNumericColumns(:, 6));

% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns idx;
%% Create output variable
SZ_tunnel.x = SZ_tunnel.x-GTS_coordinates.x;
SZ_tunnel.y = SZ_tunnel.y-GTS_coordinates.y;
SZ_tunnel.z = SZ_tunnel.z-GTS_coordinates.z;
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;
end

