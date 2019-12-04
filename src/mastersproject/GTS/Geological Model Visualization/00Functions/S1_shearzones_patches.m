function [S11, S12, S13] =S1_shearzones_patches(SZ_tunnel, FBS, SBH, PRP, INJ, GEO, iv)
%% ======================== Import Coordinates ============================
%% ======================== Import S1_1 ===================================
filename = 'S1_1.txt';
delimiter = '\t';
startRow = 2;
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
rawData = dataArray{2};
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
            numericData(row, 2) = numbers{1};
            raw{row, 2} = numbers{1};
        end
    catch
        raw{row, 2} = rawData{row};
    end
end
rawNumericColumns = raw(:, 2);
rawStringColumns = string(raw(:, 1));
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";
% Create output variable
S11 = table;
S11.Borehole = categorical(rawStringColumns(:, 1));
S11.Depthm = cell2mat(rawNumericColumns(:, 1));
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R idx;


for ii = 1 : length(S11.Borehole(:,1))
    if S11.Borehole(ii,1) == 'INJ1'
        S11.x(ii,1) = INJ(1,1)+S11.Depthm(ii).*sind(INJ(1,6)).*cosd(INJ(1,7));
        S11.y(ii,1) = INJ(1,2)+S11.Depthm(ii).*cosd(INJ(1,6))*cosd(INJ(1,7));
        S11.z(ii,1) = INJ(1,3)+S11.Depthm(ii).*sind(INJ(1,7));
    elseif S11.Borehole(ii,1) == 'INJ2'
        S11.x(ii,1) = [INJ(2,1)+S11.Depthm(ii).*sind(INJ(2,6)).*cosd(INJ(2,7))];
        S11.y(ii,1) = [INJ(2,2)+S11.Depthm(ii).*cosd(INJ(2,6))*cosd(INJ(2,7))];
        S11.z(ii,1) = [INJ(2,3)+S11.Depthm(ii).*sind(INJ(2,7))];
    elseif S11.Borehole(ii,1) == 'FBS1'
        S11.x(ii,1) = [FBS(1,1)+S11.Depthm(ii).*sind(FBS(1,6)).*cosd(FBS(1,7))];
        S11.y(ii,1) = [FBS(1,2)+S11.Depthm(ii).*cosd(FBS(1,6))*cosd(FBS(1,7))];
        S11.z(ii,1) = [FBS(1,3)+S11.Depthm(ii).*sind(FBS(1,7))];
    elseif S11.Borehole(ii,1) == 'FBS2'
        S11.x(ii,1) = [FBS(2,1)+S11.Depthm(ii).*sind(FBS(2,6)).*cosd(FBS(2,7))];
        S11.y(ii,1) = [FBS(2,2)+S11.Depthm(ii).*cosd(FBS(2,6))*cosd(FBS(2,7))];
        S11.z(ii,1) = [FBS(2,3)+S11.Depthm(ii).*sind(FBS(2,7))];
    elseif S11.Borehole(ii,1) == 'FBS3'
        S11.x(ii,1) = [FBS(3,1)+S11.Depthm(ii).*sind(FBS(3,6)).*cosd(FBS(3,7))];
        S11.y(ii,1) = [FBS(3,2)+S11.Depthm(ii).*cosd(FBS(3,6))*cosd(FBS(3,7))];
        S11.z(ii,1) = [FBS(3,3)+S11.Depthm(ii).*sind(FBS(3,7))];
    elseif S11.Borehole(ii,1) == 'PRP1'
        S11.x(ii,1) = [PRP(1,1)+S11.Depthm(ii).*sind(PRP(1,6)).*cosd(PRP(1,7))];
        S11.y(ii,1) = [PRP(1,2)+S11.Depthm(ii).*cosd(PRP(1,6))*cosd(PRP(1,7))];
        S11.z(ii,1) = [PRP(1,3)+S11.Depthm(ii).*sind(PRP(1,7))];
    elseif S11.Borehole(ii,1) == 'PRP2'
        S11.x(ii,1) = [PRP(2,1)+S11.Depthm(ii).*sind(PRP(2,6)).*cosd(PRP(2,7))];
        S11.y(ii,1) = [PRP(2,2)+S11.Depthm(ii).*cosd(PRP(2,6))*cosd(PRP(2,7))];
        S11.z(ii,1) = [PRP(2,3)+S11.Depthm(ii).*sind(PRP(2,7))];
    elseif S11.Borehole(ii,1) == 'PRP3'
        S11.x(ii,1) = [PRP(3,1)+S11.Depthm(ii).*sind(PRP(3,6)).*cosd(PRP(3,7))];
        S11.y(ii,1) = [PRP(3,2)+S11.Depthm(ii).*cosd(PRP(3,6))*cosd(PRP(3,7))];
        S11.z(ii,1) = [PRP(3,3)+S11.Depthm(ii).*sind(PRP(3,7))];
    elseif S11.Borehole(ii,1) == 'SBH1'
        S11.x(ii,1) = [SBH(1,1)+S11.Depthm(ii).*sind(SBH(1,6)).*cosd(SBH(1,7))];
        S11.y(ii,1) = [SBH(1,2)+S11.Depthm(ii).*cosd(SBH(1,6))*cosd(SBH(1,7))];
        S11.z(ii,1) = [SBH(1,3)+S11.Depthm(ii).*sind(SBH(1,7))];
    elseif S11.Borehole(ii,1) == 'SBH3'
        S11.x(ii,1) = [SBH(2,1)+S11.Depthm(ii).*sind(SBH(2,6)).*cosd(SBH(2,7))];
        S11.y(ii,1) = [SBH(2,2)+S11.Depthm(ii).*cosd(SBH(2,6))*cosd(SBH(2,7))];
        S11.z(ii,1) = [SBH(2,3)+S11.Depthm(ii).*sind(SBH(2,7))];
    elseif S11.Borehole(ii,1) == 'SBH4'
        S11.x(ii,1) = [SBH(3,1)+S11.Depthm(ii).*sind(SBH(3,6)).*cosd(SBH(3,7))];
        S11.y(ii,1) = [SBH(3,2)+S11.Depthm(ii).*cosd(SBH(3,6))*cosd(SBH(3,7))];
        S11.z(ii,1) = [SBH(3,3)+S11.Depthm(ii).*sind(SBH(3,7))];
    elseif S11.Borehole(ii,1) == 'GEO1'
        S11.x(ii,1) = [GEO(1,1)+S11.Depthm(ii).*sind(GEO(1,6)).*cosd(GEO(1,7))];
        S11.y(ii,1) = [GEO(1,2)+S11.Depthm(ii).*cosd(GEO(1,6))*cosd(GEO(1,7))];
        S11.z(ii,1) = [GEO(1,3)+S11.Depthm(ii).*sind(GEO(1,7))];
    elseif S11.Borehole(ii,1) == 'GEO2'
        S11.x(ii,1) = [GEO(2,1)+S11.Depthm(ii).*sind(GEO(2,6)).*cosd(GEO(2,7))];
        S11.y(ii,1) = [GEO(2,2)+S11.Depthm(ii).*cosd(GEO(2,6))*cosd(GEO(2,7))];
        S11.z(ii,1) = [GEO(2,3)+S11.Depthm(ii).*sind(GEO(2,7))];
    elseif S11.Borehole(ii,1) == 'GEO3'
        S11.x(ii,1) = [GEO(3,1)+S11.Depthm(ii).*sind(GEO(3,6)).*cosd(GEO(3,7))];
        S11.y(ii,1) = [GEO(3,2)+S11.Depthm(ii).*cosd(GEO(3,6))*cosd(GEO(3,7))];
        S11.z(ii,1) = [GEO(3,3)+S11.Depthm(ii).*sind(GEO(3,7))];
    elseif S11.Borehole(ii,1) == 'GEO4'
        S11.x(ii,1) = [GEO(4,1)+S11.Depthm(ii).*sind(GEO(4,6)).*cosd(GEO(4,7))];
        S11.y(ii,1) = [GEO(4,2)+S11.Depthm(ii).*cosd(GEO(4,6))*cosd(GEO(4,7))];
        S11.z(ii,1) = [GEO(4,3)+S11.Depthm(ii).*sind(GEO(4,7))];
    end
end


%% ================= S1_2 =================================================
filename = 'S1_2.txt';
delimiter = '\t';
startRow = 2;
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
rawData = dataArray{2};
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
            numericData(row, 2) = numbers{1};
            raw{row, 2} = numbers{1};
        end
    catch
        raw{row, 2} = rawData{row};
    end
end
rawNumericColumns = raw(:, 2);
rawStringColumns = string(raw(:, 1));
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";
% Create output variable
S12 = table;
S12.Borehole = categorical(rawStringColumns(:, 1));
S12.Depthm = cell2mat(rawNumericColumns(:, 1));
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R idx;


for ii = 1 : length(S12.Borehole(:,1))
    if S12.Borehole(ii,1) == 'INJ1'
        S12.x(ii,1) = INJ(1,1)+S12.Depthm(ii).*sind(INJ(1,6)).*cosd(INJ(1,7));
        S12.y(ii,1) = INJ(1,2)+S12.Depthm(ii).*cosd(INJ(1,6))*cosd(INJ(1,7));
        S12.z(ii,1) = INJ(1,3)+S12.Depthm(ii).*sind(INJ(1,7));
    elseif S12.Borehole(ii,1) == 'INJ2'
        S12.x(ii,1) = [INJ(2,1)+S12.Depthm(ii).*sind(INJ(2,6)).*cosd(INJ(2,7))];
        S12.y(ii,1) = [INJ(2,2)+S12.Depthm(ii).*cosd(INJ(2,6))*cosd(INJ(2,7))];
        S12.z(ii,1) = [INJ(2,3)+S12.Depthm(ii).*sind(INJ(2,7))];
    elseif S12.Borehole(ii,1) == 'FBS1'
        S12.x(ii,1) = [FBS(1,1)+S12.Depthm(ii).*sind(FBS(1,6)).*cosd(FBS(1,7))];
        S12.y(ii,1) = [FBS(1,2)+S12.Depthm(ii).*cosd(FBS(1,6))*cosd(FBS(1,7))];
        S12.z(ii,1) = [FBS(1,3)+S12.Depthm(ii).*sind(FBS(1,7))];
    elseif S12.Borehole(ii,1) == 'FBS2'
        S12.x(ii,1) = [FBS(2,1)+S12.Depthm(ii).*sind(FBS(2,6)).*cosd(FBS(2,7))];
        S12.y(ii,1) = [FBS(2,2)+S12.Depthm(ii).*cosd(FBS(2,6))*cosd(FBS(2,7))];
        S12.z(ii,1) = [FBS(2,3)+S12.Depthm(ii).*sind(FBS(2,7))];
    elseif S12.Borehole(ii,1) == 'FBS3'
        S12.x(ii,1) = [FBS(3,1)+S12.Depthm(ii).*sind(FBS(3,6)).*cosd(FBS(3,7))];
        S12.y(ii,1) = [FBS(3,2)+S12.Depthm(ii).*cosd(FBS(3,6))*cosd(FBS(3,7))];
        S12.z(ii,1) = [FBS(3,3)+S12.Depthm(ii).*sind(FBS(3,7))];
    elseif S12.Borehole(ii,1) == 'PRP1'
        S12.x(ii,1) = [PRP(1,1)+S12.Depthm(ii).*sind(PRP(1,6)).*cosd(PRP(1,7))];
        S12.y(ii,1) = [PRP(1,2)+S12.Depthm(ii).*cosd(PRP(1,6))*cosd(PRP(1,7))];
        S12.z(ii,1) = [PRP(1,3)+S12.Depthm(ii).*sind(PRP(1,7))];
    elseif S12.Borehole(ii,1) == 'PRP2'
        S12.x(ii,1) = [PRP(2,1)+S12.Depthm(ii).*sind(PRP(2,6)).*cosd(PRP(2,7))];
        S12.y(ii,1) = [PRP(2,2)+S12.Depthm(ii).*cosd(PRP(2,6))*cosd(PRP(2,7))];
        S12.z(ii,1) = [PRP(2,3)+S12.Depthm(ii).*sind(PRP(2,7))];
    elseif S12.Borehole(ii,1) == 'PRP3'
        S12.x(ii,1) = [PRP(3,1)+S12.Depthm(ii).*sind(PRP(3,6)).*cosd(PRP(3,7))];
        S12.y(ii,1) = [PRP(3,2)+S12.Depthm(ii).*cosd(PRP(3,6))*cosd(PRP(3,7))];
        S12.z(ii,1) = [PRP(3,3)+S12.Depthm(ii).*sind(PRP(3,7))];
    elseif S12.Borehole(ii,1) == 'SBH1'
        S12.x(ii,1) = [SBH(1,1)+S12.Depthm(ii).*sind(SBH(1,6)).*cosd(SBH(1,7))];
        S12.y(ii,1) = [SBH(1,2)+S12.Depthm(ii).*cosd(SBH(1,6))*cosd(SBH(1,7))];
        S12.z(ii,1) = [SBH(1,3)+S12.Depthm(ii).*sind(SBH(1,7))];
    elseif S12.Borehole(ii,1) == 'SBH3'
        S12.x(ii,1) = [SBH(2,1)+S12.Depthm(ii).*sind(SBH(2,6)).*cosd(SBH(2,7))];
        S12.y(ii,1) = [SBH(2,2)+S12.Depthm(ii).*cosd(SBH(2,6))*cosd(SBH(2,7))];
        S12.z(ii,1) = [SBH(2,3)+S12.Depthm(ii).*sind(SBH(2,7))];
    elseif S12.Borehole(ii,1) == 'SBH4'
        S12.x(ii,1) = [SBH(3,1)+S12.Depthm(ii).*sind(SBH(3,6)).*cosd(SBH(3,7))];
        S12.y(ii,1) = [SBH(3,2)+S12.Depthm(ii).*cosd(SBH(3,6))*cosd(SBH(3,7))];
        S12.z(ii,1) = [SBH(3,3)+S12.Depthm(ii).*sind(SBH(3,7))];
    elseif S12.Borehole(ii,1) == 'GEO1'
        S12.x(ii,1) = [GEO(1,1)+S12.Depthm(ii).*sind(GEO(1,6)).*cosd(GEO(1,7))];
        S12.y(ii,1) = [GEO(1,2)+S12.Depthm(ii).*cosd(GEO(1,6))*cosd(GEO(1,7))];
        S12.z(ii,1) = [GEO(1,3)+S12.Depthm(ii).*sind(GEO(1,7))];
    elseif S12.Borehole(ii,1) == 'GEO2'
        S12.x(ii,1) = [GEO(2,1)+S12.Depthm(ii).*sind(GEO(2,6)).*cosd(GEO(2,7))];
        S12.y(ii,1) = [GEO(2,2)+S12.Depthm(ii).*cosd(GEO(2,6))*cosd(GEO(2,7))];
        S12.z(ii,1) = [GEO(2,3)+S12.Depthm(ii).*sind(GEO(2,7))];
    elseif S12.Borehole(ii,1) == 'GEO3'
        S12.x(ii,1) = [GEO(3,1)+S12.Depthm(ii).*sind(GEO(3,6)).*cosd(GEO(3,7))];
        S12.y(ii,1) = [GEO(3,2)+S12.Depthm(ii).*cosd(GEO(3,6))*cosd(GEO(3,7))];
        S12.z(ii,1) = [GEO(3,3)+S12.Depthm(ii).*sind(GEO(3,7))];
    elseif S12.Borehole(ii,1) == 'GEO4'
        S12.x(ii,1) = [GEO(4,1)+S12.Depthm(ii).*sind(GEO(4,6)).*cosd(GEO(4,7))];
        S12.y(ii,1) = [GEO(4,2)+S12.Depthm(ii).*cosd(GEO(4,6))*cosd(GEO(4,7))];
        S12.z(ii,1) = [GEO(4,3)+S12.Depthm(ii).*sind(GEO(4,7))];
    end
end

%% ======================== S1_3 ==========================================
filename = 'S1_3.txt';
delimiter = '\t';
startRow = 2;
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
rawData = dataArray{2};
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
            numericData(row, 2) = numbers{1};
            raw{row, 2} = numbers{1};
        end
    catch
        raw{row, 2} = rawData{row};
    end
end
rawNumericColumns = raw(:, 2);
rawStringColumns = string(raw(:, 1));
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
idx = (rawStringColumns(:, 1) == "<undefined>");
rawStringColumns(idx, 1) = "";
% Create output variable
S13 = table;
S13.Borehole = categorical(rawStringColumns(:, 1));
S13.Depthm = cell2mat(rawNumericColumns(:, 1));
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R idx;


for ii = 1 : length(S13.Borehole(:,1))
    if S13.Borehole(ii,1) == 'INJ1'
        S13.x(ii,1) = INJ(1,1)+S13.Depthm(ii).*sind(INJ(1,6)).*cosd(INJ(1,7));
        S13.y(ii,1) = INJ(1,2)+S13.Depthm(ii).*cosd(INJ(1,6))*cosd(INJ(1,7));
        S13.z(ii,1) = INJ(1,3)+S13.Depthm(ii).*sind(INJ(1,7));
    elseif S13.Borehole(ii,1) == 'INJ2'
        S13.x(ii,1) = [INJ(2,1)+S13.Depthm(ii).*sind(INJ(2,6)).*cosd(INJ(2,7))];
        S13.y(ii,1) = [INJ(2,2)+S13.Depthm(ii).*cosd(INJ(2,6))*cosd(INJ(2,7))];
        S13.z(ii,1) = [INJ(2,3)+S13.Depthm(ii).*sind(INJ(2,7))];
    elseif S13.Borehole(ii,1) == 'FBS1'
        S13.x(ii,1) = [FBS(1,1)+S13.Depthm(ii).*sind(FBS(1,6)).*cosd(FBS(1,7))];
        S13.y(ii,1) = [FBS(1,2)+S13.Depthm(ii).*cosd(FBS(1,6))*cosd(FBS(1,7))];
        S13.z(ii,1) = [FBS(1,3)+S13.Depthm(ii).*sind(FBS(1,7))];
    elseif S13.Borehole(ii,1) == 'FBS2'
        S13.x(ii,1) = [FBS(2,1)+S13.Depthm(ii).*sind(FBS(2,6)).*cosd(FBS(2,7))];
        S13.y(ii,1) = [FBS(2,2)+S13.Depthm(ii).*cosd(FBS(2,6))*cosd(FBS(2,7))];
        S13.z(ii,1) = [FBS(2,3)+S13.Depthm(ii).*sind(FBS(2,7))];
    elseif S13.Borehole(ii,1) == 'FBS3'
        S13.x(ii,1) = [FBS(3,1)+S13.Depthm(ii).*sind(FBS(3,6)).*cosd(FBS(3,7))];
        S13.y(ii,1) = [FBS(3,2)+S13.Depthm(ii).*cosd(FBS(3,6))*cosd(FBS(3,7))];
        S13.z(ii,1) = [FBS(3,3)+S13.Depthm(ii).*sind(FBS(3,7))];
    elseif S13.Borehole(ii,1) == 'PRP1'
        S13.x(ii,1) = [PRP(1,1)+S13.Depthm(ii).*sind(PRP(1,6)).*cosd(PRP(1,7))];
        S13.y(ii,1) = [PRP(1,2)+S13.Depthm(ii).*cosd(PRP(1,6))*cosd(PRP(1,7))];
        S13.z(ii,1) = [PRP(1,3)+S13.Depthm(ii).*sind(PRP(1,7))];
    elseif S13.Borehole(ii,1) == 'PRP2'
        S13.x(ii,1) = [PRP(2,1)+S13.Depthm(ii).*sind(PRP(2,6)).*cosd(PRP(2,7))];
        S13.y(ii,1) = [PRP(2,2)+S13.Depthm(ii).*cosd(PRP(2,6))*cosd(PRP(2,7))];
        S13.z(ii,1) = [PRP(2,3)+S13.Depthm(ii).*sind(PRP(2,7))];
    elseif S13.Borehole(ii,1) == 'PRP3'
        S13.x(ii,1) = [PRP(3,1)+S13.Depthm(ii).*sind(PRP(3,6)).*cosd(PRP(3,7))];
        S13.y(ii,1) = [PRP(3,2)+S13.Depthm(ii).*cosd(PRP(3,6))*cosd(PRP(3,7))];
        S13.z(ii,1) = [PRP(3,3)+S13.Depthm(ii).*sind(PRP(3,7))];
    elseif S13.Borehole(ii,1) == 'SBH1'
        S13.x(ii,1) = [SBH(1,1)+S13.Depthm(ii).*sind(SBH(1,6)).*cosd(SBH(1,7))];
        S13.y(ii,1) = [SBH(1,2)+S13.Depthm(ii).*cosd(SBH(1,6))*cosd(SBH(1,7))];
        S13.z(ii,1) = [SBH(1,3)+S13.Depthm(ii).*sind(SBH(1,7))];
    elseif S13.Borehole(ii,1) == 'SBH3'
        S13.x(ii,1) = [SBH(2,1)+S13.Depthm(ii).*sind(SBH(2,6)).*cosd(SBH(2,7))];
        S13.y(ii,1) = [SBH(2,2)+S13.Depthm(ii).*cosd(SBH(2,6))*cosd(SBH(2,7))];
        S13.z(ii,1) = [SBH(2,3)+S13.Depthm(ii).*sind(SBH(2,7))];
    elseif S13.Borehole(ii,1) == 'SBH4'
        S13.x(ii,1) = [SBH(3,1)+S13.Depthm(ii).*sind(SBH(3,6)).*cosd(SBH(3,7))];
        S13.y(ii,1) = [SBH(3,2)+S13.Depthm(ii).*cosd(SBH(3,6))*cosd(SBH(3,7))];
        S13.z(ii,1) = [SBH(3,3)+S13.Depthm(ii).*sind(SBH(3,7))];
    elseif S13.Borehole(ii,1) == 'GEO1'
        S13.x(ii,1) = [GEO(1,1)+S13.Depthm(ii).*sind(GEO(1,6)).*cosd(GEO(1,7))];
        S13.y(ii,1) = [GEO(1,2)+S13.Depthm(ii).*cosd(GEO(1,6))*cosd(GEO(1,7))];
        S13.z(ii,1) = [GEO(1,3)+S13.Depthm(ii).*sind(GEO(1,7))];
    elseif S13.Borehole(ii,1) == 'GEO2'
        S13.x(ii,1) = [GEO(2,1)+S13.Depthm(ii).*sind(GEO(2,6)).*cosd(GEO(2,7))];
        S13.y(ii,1) = [GEO(2,2)+S13.Depthm(ii).*cosd(GEO(2,6))*cosd(GEO(2,7))];
        S13.z(ii,1) = [GEO(2,3)+S13.Depthm(ii).*sind(GEO(2,7))];
    elseif S13.Borehole(ii,1) == 'GEO3'
        S13.x(ii,1) = [GEO(3,1)+S13.Depthm(ii).*sind(GEO(3,6)).*cosd(GEO(3,7))];
        S13.y(ii,1) = [GEO(3,2)+S13.Depthm(ii).*cosd(GEO(3,6))*cosd(GEO(3,7))];
        S13.z(ii,1) = [GEO(3,3)+S13.Depthm(ii).*sind(GEO(3,7))];
    elseif S13.Borehole(ii,1) == 'GEO4'
        S13.x(ii,1) = [GEO(4,1)+S13.Depthm(ii).*sind(GEO(4,6)).*cosd(GEO(4,7))];
        S13.y(ii,1) = [GEO(4,2)+S13.Depthm(ii).*cosd(GEO(4,6))*cosd(GEO(4,7))];
        S13.z(ii,1) = [GEO(4,3)+S13.Depthm(ii).*sind(GEO(4,7))];
    end
end

%% =============== Assign tunnel intersections to Shear zones =============
for ii = 1:length(SZ_tunnel.Tunnel)
if SZ_tunnel.SZ(ii) == 11
    S11.Borehole(end+1,1) = SZ_tunnel.Tunnel(ii);
    S11.Depthm(end,1) = NaN;
    S11.x(end,1) = SZ_tunnel.x(ii);
    S11.y(end,1) = SZ_tunnel.y(ii);
    S11.z(end,1) = SZ_tunnel.z(ii);
elseif SZ_tunnel.SZ(ii) == 12
        S12.Borehole(end+1,1) = SZ_tunnel.Tunnel(ii);
    S12.Depthm(end,1) = NaN;
    S12.x(end,1) = SZ_tunnel.x(ii);
    S12.y(end,1) = SZ_tunnel.y(ii);
    S12.z(end,1) = SZ_tunnel.z(ii);
elseif SZ_tunnel.SZ(ii) == 13
        S13.Borehole(end+1,1) = SZ_tunnel.Tunnel(ii);
    S13.Depthm(end,1) = NaN;
    S13.x(end,1) = SZ_tunnel.x(ii);
    S13.y(end,1) = SZ_tunnel.y(ii);
    S13.z(end,1) = SZ_tunnel.z(ii);
end
end


%% ============================ Build patches =============================

figure(iv)
% S1_1 Patches
patches.S11(1,1:3) = [{'VE'}, {'FBS1'}, {'INJ1'}];
patches.S11(2,1:3) = [{'INJ1'},{'FBS1'},{'PRP2'}];
patches.S11(3,1:3) = [{'INJ1'},{'FBS1'},{'PRP1'}];
patches.S11(4,1:3) = [{'PRP1'},{'PRP2'},{'INJ2'}];
patches.S11(5,1:3) = [{'AU'},{'PRP3'},{'INJ2'}];
patches.S11(6,1:3) = [{'PRP3'},{'INJ2'},{'PRP2'}];
patches.S11(7,1:3) = [{'PRP2'},{'FBS1'},{'PRP1'}];
patches.S11(8,1:3) = [{'PRP3'},{'GEO3'},{'FBS2'}];
patches.S11(9,1:3) = [{'PRP3'},{'GEO3'},{'GEO4'}];
patches.S11(10,1:3) = [{'PRP2'},{'GEO3'},{'GEO4'}];
patches.S11(11,1:3) = [{'PRP3'},{'AU'},{'GEO4'}]; 
patches.S11(12,1:3) = [{'FBS1'},{'VE'},{'PRP2'}]; 
patches.S11(13,1:3) = [{'GEO4'},{'VE'},{'FBS1'}]; 
patches.S11(14,1:3) = [{'GEO4'},{'VE'},{'AU'}]; 

for ii = 1: length(patches.S11(:,1)) 
    for iii = 1:length(patches.S11(ii,:))
       loc = find(patches.S11(ii,iii) == S11.Borehole);
        x_patch(ii,iii) = S11.x(loc);
        y_patch(ii,iii) = S11.y(loc);
        z_patch(ii,iii) = S11.z(loc);
    end
end
for ii = 1:length(x_patch(:,1))
    X = x_patch(ii,:);
    Y = y_patch(ii,:);
    Z = z_patch(ii,:);
    patch(X,Y,Z,'r','LineStyle','none','FaceAlpha',0.7);
end

clear X Y Z x_patch y_patch z_patch ii iii
% S1_2 Patches
patches.S12(1,1:3) = [{'AU'}, {'PRP3'}, {'INJ2'}];
patches.S12(2,1:3) = [{'INJ2'},{'PRP3'},{'GEO3'}];
patches.S12(3,1:3) = [{'INJ2'},{'GEO3'},{'GEO4'}];
patches.S12(4,1:3) = [{'INJ2'},{'INJ1'},{'GEO4'}];
patches.S12(5,1:3) = [{'VE'},{'GEO4'},{'INJ1'}];
patches.S12(6,1:3) = [{'PRP3'},{'INJ2'},{'GEO4'}];
patches.S12(7,1:3) = [{'FBS1'},{'INJ1'},{'GEO4'}]; 
patches.S12(8,1:3) = [{'FBS1'},{'INJ2'},{'GEO4'}]; 
patches.S12(9,1:3) = [{'PRP3'},{'AU'},{'GEO4'}]; 
patches.S12(10,1:3) = [{'GEO4'},{'VE'},{'AU'}]; 


for ii = 1: length(patches.S12(:,1)) 
    for iii = 1:length(patches.S12(ii,:))
       loc = find(patches.S12(ii,iii) == S12.Borehole);
        x_patch(ii,iii) = S12.x(loc);
        y_patch(ii,iii) = S12.y(loc);
        z_patch(ii,iii) = S12.z(loc);
    end
end
for ii = 1:length(x_patch(:,1))
    X = x_patch(ii,:);
    Y = y_patch(ii,:);
    Z = z_patch(ii,:);
    patch(X,Y,Z,'r','LineStyle','none','FaceAlpha',0.7);
end


clear X Y Z x_patch y_patch z_patch ii iii
%S1_3 Patches
% patches.S13(1,1:3) = [{'AU'}, {'PRP1'}, {'PRP2'}];
patches.S13(1,1:3) = [{'PRP1'},{'INJ2'},{'FBS1'}];
patches.S13(2,1:3) = [{'FBS1'},{'PRP1'},{'VE'}];
patches.S13(3,1:3) = [{'PRP2'},{'FBS1'},{'VE'}];
patches.S13(4,1:3) = [{'AU'},{'PRP2'},{'INJ2'}];
patches.S13(5,1:3) = [{'FBS1'},{'INJ2'},{'PRP2'}];
patches.S13(6,1:3) = [{'AU'},{'VE'},{'PRP2'}];         


for ii = 1: length(patches.S13(:,1)) 
    for iii = 1:length(patches.S13(ii,:))
       loc = find(patches.S13(ii,iii) == S13.Borehole);
        x_patch(ii,iii) = S13.x(loc);
        y_patch(ii,iii) = S13.y(loc);
        z_patch(ii,iii) = S13.z(loc);
    end
end
for ii = 1:length(x_patch(:,1))
    X = x_patch(ii,:);
    Y = y_patch(ii,:);
    Z = z_patch(ii,:);
    patch(X,Y,Z,'r','LineStyle','none','FaceAlpha',0.7);
end

















end