function [SBH1, SBH3, SBH4] = Stress_measurement_locations(SBH_coord)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


filename = '..\09_Monitoring\02_stress_measurement_locations\Stress_Locations.txt';
delimiter = '\t';

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2]
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
rawNumericColumns = raw(:, [1,2]);
rawStringColumns = string(raw(:, 3));
%% Create output variable
Locations = table;
Locations.Depth = cell2mat(rawNumericColumns(:, 1));
Locations.Borehole = cell2mat(rawNumericColumns(:, 2));
Locations.Method = rawStringColumns(:, 1);

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns;
for iip = 3:4

count1 = 1;
count3 = 3;
count4 = 4;
for ii = 1:length(Locations.Depth(:,1))
    if Locations.Borehole(ii)== 1
        SBH1.loc(count1,1:3) = [SBH_coord(1,1)+Locations.Depth(ii).*sind(258.89).*cosd(-75.13), SBH_coord(1,2)+Locations.Depth(ii).*cosd(258.89)*cosd(-75.13), SBH_coord(1,3)+Locations.Depth(ii).*sind(-75.13)];
        SBH1.meth(count1,:) = Locations.Method(ii);
        count1 = count1+1;
    elseif Locations.Borehole(ii)== 3
        SBH3.loc(count3,1:3) = [SBH_coord(3,1)+Locations.Depth(ii).*sind(190.07).*cosd(4.91), SBH_coord(3,2)+Locations.Depth(ii).*cosd(190.07)*cosd(4.91), SBH_coord(3,3)+Locations.Depth(ii).*sind(4.91)];
        SBH3.meth(count3,:) = Locations.Method(ii);
        count3 = count3+1;
    elseif Locations.Borehole(ii)== 4
        SBH4.loc(count4,1:3) = [SBH_coord(5,1)+Locations.Depth(ii).*sind(320).*cosd(5), SBH_coord(5,2)+Locations.Depth(ii).*cosd(320)*cosd(5), SBH_coord(5,3)+Locations.Depth(ii).*sind(5)];
        SBH4.meth(count4,:) = Locations.Method(ii);
        count4 = count4+1;
    end
end

for ii = 1:length(SBH1.loc(:,1))
    if strcmp(SBH1.meth(ii),{'USBM'})
        color = 'r';
    elseif strcmp(SBH1.meth(ii),{'CSIRO'})
        color = 'g';
    elseif strcmp(SBH1.meth(ii),{'HF'})
        color = 'b';
    end

    scatter3(SBH1.loc(ii,1),SBH1.loc(ii,2),SBH1.loc(ii,3),[],color,'filled')
end


for ii = 1:length(SBH3.loc(:,1))
    if strcmp(SBH3.meth(ii),{'USBM'})
        color = 'r';
    elseif strcmp(SBH3.meth(ii),{'CSIRO'})
        color = 'g';
    elseif strcmp(SBH3.meth(ii),{'HF'})
        color = 'b';
    end

    scatter3(SBH3.loc(ii,1),SBH3.loc(ii,2),SBH3.loc(ii,3),[],color,'filled')
end

for ii = 1:length(SBH4.loc(:,1))
    if strcmp(SBH4.meth(ii),{'USBM'})
        color = 'r';
    elseif strcmp(SBH4.meth(ii),{'CSIRO'})
        color = 'g';
    elseif strcmp(SBH4.meth(ii),{'HF'})
        color = 'b';
    end

    scatter3(SBH4.loc(ii,1),SBH4.loc(ii,2),SBH4.loc(ii,3),[],color,'filled')
end
end
end
