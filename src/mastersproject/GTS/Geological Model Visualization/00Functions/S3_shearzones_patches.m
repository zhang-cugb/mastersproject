function [S31, S32] =S3_shearzones_patches(SZ_tunnel, FBS, SBH, PRP, INJ, GEO,iv)
%% ======================== Import Coordinates ============================
%% ======================== Import S3_1 ===================================
filename = 'S3_1.txt';
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
S31 = table;
S31.Borehole = categorical(rawStringColumns(:, 1));
S31.Depthm = cell2mat(rawNumericColumns(:, 1));
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R idx;


for ii = 1 : length(S31.Borehole(:,1))
   if S31.Borehole(ii,1) == 'INJ1'
        S31.x(ii,1) = INJ(1,1)+S31.Depthm(ii).*sind(INJ(1,6)).*cosd(INJ(1,7));
        S31.y(ii,1) = INJ(1,2)+S31.Depthm(ii).*cosd(INJ(1,6))*cosd(INJ(1,7));
        S31.z(ii,1) = INJ(1,3)+S31.Depthm(ii).*sind(INJ(1,7));
    elseif S31.Borehole(ii,1) == 'INJ2'
        S31.x(ii,1) = [INJ(2,1)+S31.Depthm(ii).*sind(INJ(2,6)).*cosd(INJ(2,7))];
        S31.y(ii,1) = [INJ(2,2)+S31.Depthm(ii).*cosd(INJ(2,6))*cosd(INJ(2,7))];
        S31.z(ii,1) = [INJ(2,3)+S31.Depthm(ii).*sind(INJ(2,7))];
    elseif S31.Borehole(ii,1) == 'FBS1'
        S31.x(ii,1) = [FBS(1,1)+S31.Depthm(ii).*sind(FBS(1,6)).*cosd(FBS(1,7))];
        S31.y(ii,1) = [FBS(1,2)+S31.Depthm(ii).*cosd(FBS(1,6))*cosd(FBS(1,7))];
        S31.z(ii,1) = [FBS(1,3)+S31.Depthm(ii).*sind(FBS(1,7))];
    elseif S31.Borehole(ii,1) == 'FBS2'
        S31.x(ii,1) = [FBS(2,1)+S31.Depthm(ii).*sind(FBS(2,6)).*cosd(FBS(2,7))];
        S31.y(ii,1) = [FBS(2,2)+S31.Depthm(ii).*cosd(FBS(2,6))*cosd(FBS(2,7))];
        S31.z(ii,1) = [FBS(2,3)+S31.Depthm(ii).*sind(FBS(2,7))];
    elseif S31.Borehole(ii,1) == 'FBS3'
        S31.x(ii,1) = [FBS(3,1)+S31.Depthm(ii).*sind(FBS(3,6)).*cosd(FBS(3,7))];
        S31.y(ii,1) = [FBS(3,2)+S31.Depthm(ii).*cosd(FBS(3,6))*cosd(FBS(3,7))];
        S31.z(ii,1) = [FBS(3,3)+S31.Depthm(ii).*sind(FBS(3,7))];
    elseif S31.Borehole(ii,1) == 'PRP1'
        S31.x(ii,1) = [PRP(1,1)+S31.Depthm(ii).*sind(PRP(1,6)).*cosd(PRP(1,7))];
        S31.y(ii,1) = [PRP(1,2)+S31.Depthm(ii).*cosd(PRP(1,6))*cosd(PRP(1,7))];
        S31.z(ii,1) = [PRP(1,3)+S31.Depthm(ii).*sind(PRP(1,7))];
    elseif S31.Borehole(ii,1) == 'PRP2'
        S31.x(ii,1) = [PRP(2,1)+S31.Depthm(ii).*sind(PRP(2,6)).*cosd(PRP(2,7))];
        S31.y(ii,1) = [PRP(2,2)+S31.Depthm(ii).*cosd(PRP(2,6))*cosd(PRP(2,7))];
        S31.z(ii,1) = [PRP(2,3)+S31.Depthm(ii).*sind(PRP(2,7))];
    elseif S31.Borehole(ii,1) == 'PRP3'
        S31.x(ii,1) = [PRP(3,1)+S31.Depthm(ii).*sind(PRP(3,6)).*cosd(PRP(3,7))];
        S31.y(ii,1) = [PRP(3,2)+S31.Depthm(ii).*cosd(PRP(3,6))*cosd(PRP(3,7))];
        S31.z(ii,1) = [PRP(3,3)+S31.Depthm(ii).*sind(PRP(3,7))];
    elseif S31.Borehole(ii,1) == 'SBH1'
        S31.x(ii,1) = [SBH(1,1)+S31.Depthm(ii).*sind(SBH(1,6)).*cosd(SBH(1,7))];
        S31.y(ii,1) = [SBH(1,2)+S31.Depthm(ii).*cosd(SBH(1,6))*cosd(SBH(1,7))];
        S31.z(ii,1) = [SBH(1,3)+S31.Depthm(ii).*sind(SBH(1,7))];
    elseif S31.Borehole(ii,1) == 'SBH3'
        S31.x(ii,1) = [SBH(2,1)+S31.Depthm(ii).*sind(SBH(2,6)).*cosd(SBH(2,7))];
        S31.y(ii,1) = [SBH(2,2)+S31.Depthm(ii).*cosd(SBH(2,6))*cosd(SBH(2,7))];
        S31.z(ii,1) = [SBH(2,3)+S31.Depthm(ii).*sind(SBH(2,7))];
    elseif S31.Borehole(ii,1) == 'SBH4'
        S31.x(ii,1) = [SBH(3,1)+S31.Depthm(ii).*sind(SBH(3,6)).*cosd(SBH(3,7))];
        S31.y(ii,1) = [SBH(3,2)+S31.Depthm(ii).*cosd(SBH(3,6))*cosd(SBH(3,7))];
        S31.z(ii,1) = [SBH(3,3)+S31.Depthm(ii).*sind(SBH(3,7))];
    elseif S31.Borehole(ii,1) == 'GEO1'
        S31.x(ii,1) = [GEO(1,1)+S31.Depthm(ii).*sind(GEO(1,6)).*cosd(GEO(1,7))];
        S31.y(ii,1) = [GEO(1,2)+S31.Depthm(ii).*cosd(GEO(1,6))*cosd(GEO(1,7))];
        S31.z(ii,1) = [GEO(1,3)+S31.Depthm(ii).*sind(GEO(1,7))];
    elseif S31.Borehole(ii,1) == 'GEO2'
        S31.x(ii,1) = [GEO(2,1)+S31.Depthm(ii).*sind(GEO(2,6)).*cosd(GEO(2,7))];
        S31.y(ii,1) = [GEO(2,2)+S31.Depthm(ii).*cosd(GEO(2,6))*cosd(GEO(2,7))];
        S31.z(ii,1) = [GEO(2,3)+S31.Depthm(ii).*sind(GEO(2,7))];
    elseif S31.Borehole(ii,1) == 'GEO3'
        S31.x(ii,1) = [GEO(3,1)+S31.Depthm(ii).*sind(GEO(3,6)).*cosd(GEO(3,7))];
        S31.y(ii,1) = [GEO(3,2)+S31.Depthm(ii).*cosd(GEO(3,6))*cosd(GEO(3,7))];
        S31.z(ii,1) = [GEO(3,3)+S31.Depthm(ii).*sind(GEO(3,7))];
    elseif S31.Borehole(ii,1) == 'GEO4'
        S31.x(ii,1) = [GEO(4,1)+S31.Depthm(ii).*sind(GEO(4,6)).*cosd(GEO(4,7))];
        S31.y(ii,1) = [GEO(4,2)+S31.Depthm(ii).*cosd(GEO(4,6))*cosd(GEO(4,7))];
        S31.z(ii,1) = [GEO(4,3)+S31.Depthm(ii).*sind(GEO(4,7))];
    end
end


%% ================= S3_2 =================================================
filename = 'S3_2.txt';
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
S32 = table;
S32.Borehole = categorical(rawStringColumns(:, 1));
S32.Depthm = cell2mat(rawNumericColumns(:, 1));
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R idx;


for ii = 1 : length(S32.Borehole(:,1))
    if S32.Borehole(ii,1) == 'INJ1'
        S32.x(ii,1) = INJ(1,1)+S32.Depthm(ii).*sind(INJ(1,6)).*cosd(INJ(1,7));
        S32.y(ii,1) = INJ(1,2)+S32.Depthm(ii).*cosd(INJ(1,6))*cosd(INJ(1,7));
        S32.z(ii,1) = INJ(1,3)+S32.Depthm(ii).*sind(INJ(1,7));
    elseif S32.Borehole(ii,1) == 'INJ2'
        S32.x(ii,1) = [INJ(2,1)+S32.Depthm(ii).*sind(INJ(2,6)).*cosd(INJ(2,7))];
        S32.y(ii,1) = [INJ(2,2)+S32.Depthm(ii).*cosd(INJ(2,6))*cosd(INJ(2,7))];
        S32.z(ii,1) = [INJ(2,3)+S32.Depthm(ii).*sind(INJ(2,7))];
    elseif S32.Borehole(ii,1) == 'FBS1'
        S32.x(ii,1) = [FBS(1,1)+S32.Depthm(ii).*sind(FBS(1,6)).*cosd(FBS(1,7))];
        S32.y(ii,1) = [FBS(1,2)+S32.Depthm(ii).*cosd(FBS(1,6))*cosd(FBS(1,7))];
        S32.z(ii,1) = [FBS(1,3)+S32.Depthm(ii).*sind(FBS(1,7))];
    elseif S32.Borehole(ii,1) == 'FBS2'
        S32.x(ii,1) = [FBS(2,1)+S32.Depthm(ii).*sind(FBS(2,6)).*cosd(FBS(2,7))];
        S32.y(ii,1) = [FBS(2,2)+S32.Depthm(ii).*cosd(FBS(2,6))*cosd(FBS(2,7))];
        S32.z(ii,1) = [FBS(2,3)+S32.Depthm(ii).*sind(FBS(2,7))];
    elseif S32.Borehole(ii,1) == 'FBS3'
        S32.x(ii,1) = [FBS(3,1)+S32.Depthm(ii).*sind(FBS(3,6)).*cosd(FBS(3,7))];
        S32.y(ii,1) = [FBS(3,2)+S32.Depthm(ii).*cosd(FBS(3,6))*cosd(FBS(3,7))];
        S32.z(ii,1) = [FBS(3,3)+S32.Depthm(ii).*sind(FBS(3,7))];
    elseif S32.Borehole(ii,1) == 'PRP1'
        S32.x(ii,1) = [PRP(1,1)+S32.Depthm(ii).*sind(PRP(1,6)).*cosd(PRP(1,7))];
        S32.y(ii,1) = [PRP(1,2)+S32.Depthm(ii).*cosd(PRP(1,6))*cosd(PRP(1,7))];
        S32.z(ii,1) = [PRP(1,3)+S32.Depthm(ii).*sind(PRP(1,7))];
    elseif S32.Borehole(ii,1) == 'PRP2'
        S32.x(ii,1) = [PRP(2,1)+S32.Depthm(ii).*sind(PRP(2,6)).*cosd(PRP(2,7))];
        S32.y(ii,1) = [PRP(2,2)+S32.Depthm(ii).*cosd(PRP(2,6))*cosd(PRP(2,7))];
        S32.z(ii,1) = [PRP(2,3)+S32.Depthm(ii).*sind(PRP(2,7))];
    elseif S32.Borehole(ii,1) == 'PRP3'
        S32.x(ii,1) = [PRP(3,1)+S32.Depthm(ii).*sind(PRP(3,6)).*cosd(PRP(3,7))];
        S32.y(ii,1) = [PRP(3,2)+S32.Depthm(ii).*cosd(PRP(3,6))*cosd(PRP(3,7))];
        S32.z(ii,1) = [PRP(3,3)+S32.Depthm(ii).*sind(PRP(3,7))];
    elseif S32.Borehole(ii,1) == 'SBH1'
        S32.x(ii,1) = [SBH(1,1)+S32.Depthm(ii).*sind(SBH(1,6)).*cosd(SBH(1,7))];
        S32.y(ii,1) = [SBH(1,2)+S32.Depthm(ii).*cosd(SBH(1,6))*cosd(SBH(1,7))];
        S32.z(ii,1) = [SBH(1,3)+S32.Depthm(ii).*sind(SBH(1,7))];
    elseif S32.Borehole(ii,1) == 'SBH3'
        S32.x(ii,1) = [SBH(2,1)+S32.Depthm(ii).*sind(SBH(2,6)).*cosd(SBH(2,7))];
        S32.y(ii,1) = [SBH(2,2)+S32.Depthm(ii).*cosd(SBH(2,6))*cosd(SBH(2,7))];
        S32.z(ii,1) = [SBH(2,3)+S32.Depthm(ii).*sind(SBH(2,7))];
    elseif S32.Borehole(ii,1) == 'SBH4'
        S32.x(ii,1) = [SBH(3,1)+S32.Depthm(ii).*sind(SBH(3,6)).*cosd(SBH(3,7))];
        S32.y(ii,1) = [SBH(3,2)+S32.Depthm(ii).*cosd(SBH(3,6))*cosd(SBH(3,7))];
        S32.z(ii,1) = [SBH(3,3)+S32.Depthm(ii).*sind(SBH(3,7))];
    elseif S32.Borehole(ii,1) == 'GEO1'
        S32.x(ii,1) = [GEO(1,1)+S32.Depthm(ii).*sind(GEO(1,6)).*cosd(GEO(1,7))];
        S32.y(ii,1) = [GEO(1,2)+S32.Depthm(ii).*cosd(GEO(1,6))*cosd(GEO(1,7))];
        S32.z(ii,1) = [GEO(1,3)+S32.Depthm(ii).*sind(GEO(1,7))];
    elseif S32.Borehole(ii,1) == 'GEO2'
        S32.x(ii,1) = [GEO(2,1)+S32.Depthm(ii).*sind(GEO(2,6)).*cosd(GEO(2,7))];
        S32.y(ii,1) = [GEO(2,2)+S32.Depthm(ii).*cosd(GEO(2,6))*cosd(GEO(2,7))];
        S32.z(ii,1) = [GEO(2,3)+S32.Depthm(ii).*sind(GEO(2,7))];
    elseif S32.Borehole(ii,1) == 'GEO3'
        S32.x(ii,1) = [GEO(3,1)+S32.Depthm(ii).*sind(GEO(3,6)).*cosd(GEO(3,7))];
        S32.y(ii,1) = [GEO(3,2)+S32.Depthm(ii).*cosd(GEO(3,6))*cosd(GEO(3,7))];
        S32.z(ii,1) = [GEO(3,3)+S32.Depthm(ii).*sind(GEO(3,7))];
    elseif S32.Borehole(ii,1) == 'GEO4'
        S32.x(ii,1) = [GEO(4,1)+S32.Depthm(ii).*sind(GEO(4,6)).*cosd(GEO(4,7))];
        S32.y(ii,1) = [GEO(4,2)+S32.Depthm(ii).*cosd(GEO(4,6))*cosd(GEO(4,7))];
        S32.z(ii,1) = [GEO(4,3)+S32.Depthm(ii).*sind(GEO(4,7))];
    end
end


%% =============== Assign tunnel intersections to Shear zones =============
for ii = 1:length(SZ_tunnel.Tunnel)
if SZ_tunnel.SZ(ii) == 31
    S31.Borehole(end+1,1) = SZ_tunnel.Tunnel(ii);
    S31.Depthm(end,1) = NaN;
    S31.x(end,1) = SZ_tunnel.x(ii);
    S31.y(end,1) = SZ_tunnel.y(ii);
    S31.z(end,1) = SZ_tunnel.z(ii);
elseif SZ_tunnel.SZ(ii) == 32
        S32.Borehole(end+1,1) = SZ_tunnel.Tunnel(ii);
    S32.Depthm(end,1) = NaN;
    S32.x(end,1) = SZ_tunnel.x(ii);
    S32.y(end,1) = SZ_tunnel.y(ii);
    S32.z(end,1) = SZ_tunnel.z(ii);
end
end


%% ============================ Build patches =============================

figure(iv)
% S3_1 Patches
patches.S31(1,1:3) = [{'VE'}, {'FBS1'}, {'INJ1'}];
patches.S31(2,1:3) = [{'INJ1'},{'FBS1'},{'PRP2'}];
patches.S31(3,1:3) = [{'INJ1'},{'FBS1'},{'PRP1'}];
patches.S31(4,1:3) = [{'PRP1'},{'PRP2'},{'INJ2'}];
patches.S31(5,1:3) = [{'AU'},{'PRP3'},{'INJ2'}];
patches.S31(6,1:3) = [{'PRP3'},{'INJ2'},{'PRP2'}];
patches.S31(7,1:3) = [{'PRP2'},{'FBS1'},{'PRP1'}];
patches.S31(8,1:3) = [{'FBS1'},{'VE'},{'PRP2'}]; 
patches.S31(9,1:3) = [{'INJ1'},{'FBS3'},{'PRP1'}];
patches.S31(10,1:3) = [{'INJ1'},{'FBS3'},{'VE'}];
patches.S31(11,1:3) = [{'INJ2'},{'FBS3'},{'PRP2'}];
patches.S31(12,1:3) = [{'AU'},{'FBS3'},{'PRP2'}];
patches.S31(13,1:3) = [{'VE'},{'PRP2'},{'PRP3'}];
patches.S31(14,1:3) = [{'VE'},{'SBH4'},{'PRP3'}];
patches.S31(15,1:3) = [{'AU'},{'SBH4'},{'PRP3'}];
for ii = 1: length(patches.S31(:,1)) 
    for iii = 1:length(patches.S31(ii,:))
       loc = find(patches.S31(ii,iii) == S31.Borehole);
        x_patch(ii,iii) = S31.x(loc);
        y_patch(ii,iii) = S31.y(loc);
        z_patch(ii,iii) = S31.z(loc);
    end
end
for ii = 1:length(x_patch(:,1))
    X = x_patch(ii,:);
    Y = y_patch(ii,:);
    Z = z_patch(ii,:);
    patch(X,Y,Z,'g','LineStyle','none','FaceAlpha',0.7);
end

clear X Y Z x_patch y_patch z_patch ii iii
% % S3_2 Patches
patches.S32(1,1:3) = [{'VE'}, {'FBS1'}, {'INJ1'}];
patches.S32(2,1:3) = [{'INJ1'},{'FBS1'},{'PRP2'}];
patches.S32(3,1:3) = [{'INJ1'},{'FBS1'},{'PRP1'}];
patches.S32(4,1:3) = [{'PRP1'},{'PRP2'},{'INJ2'}];
patches.S32(5,1:3) = [{'AU'},{'PRP3'},{'INJ2'}];
patches.S32(6,1:3) = [{'PRP3'},{'INJ2'},{'PRP2'}];
patches.S32(7,1:3) = [{'PRP2'},{'FBS1'},{'PRP1'}];
patches.S32(8,1:3) = [{'FBS1'},{'VE'},{'PRP2'}]; 
patches.S32(9,1:3) = [{'INJ1'},{'FBS3'},{'PRP1'}];
patches.S32(10,1:3) = [{'INJ1'},{'FBS3'},{'VE'}];
patches.S32(11,1:3) = [{'INJ2'},{'FBS3'},{'PRP2'}];
patches.S32(12,1:3) = [{'AU'},{'FBS3'},{'PRP2'}];
patches.S32(13,1:3) = [{'VE'},{'PRP2'},{'PRP3'}];
patches.S32(14,1:3) = [{'VE'},{'SBH4'},{'PRP3'}];
patches.S32(15,1:3) = [{'AU'},{'SBH4'},{'PRP3'}];


for ii = 1: length(patches.S32(:,1)) 
    for iii = 1:length(patches.S32(ii,:))
       loc = find(patches.S32(ii,iii) == S32.Borehole);
        x_patch(ii,iii) = S32.x(loc);
        y_patch(ii,iii) = S32.y(loc);
        z_patch(ii,iii) = S32.z(loc);
    end
end
for ii = 1:length(x_patch(:,1))
    X = x_patch(ii,:);
    Y = y_patch(ii,:);
    Z = z_patch(ii,:);
    patch(X,Y,Z,'g','LineStyle','none','FaceAlpha',0.7);
end


clear X Y Z x_patch y_patch z_patch ii iii
% end

















end