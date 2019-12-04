function [ OPTVLogs ] = importOPTVstructures
% coordinates represent the midpoint of the structures.

% Choose Boreholes to import data.
boreholenames =  [{'FBS1'},{'FBS2'},{'FBS3'},{'SBH1'},{'SBH3'},{'SBH4'},...
                {'INJ1'},{'INJ2'},{'PRP1'},{'PRP2'},{'PRP3'}...
                 {'GEO1'},{'GEO2'},{'GEO3'},{'GEO4'}];
             
             
for ii = 1:length(boreholenames)
boreholename = boreholenames{ii};

[structures] = importOPTVlog( boreholename); 
OPTVLogs(1).(boreholename) = structures;

% Filter for geological structures
a = table2array(structures(:,5));
x1 = find(strcmp('Fracture',a));
x2 = find(strcmp('S1 Shear-zone',a));
x3 = find(strcmp('S3 Shear-zone',a));

OPTVLogs(2).(boreholename) = table2array(structures(x1,1:3));  %Fractures
OPTVLogs(3).(boreholename) = table2array(structures(x2,1:3));  %S1 Shear zones
OPTVLogs(4).(boreholename) = table2array(structures(x3,1:3));  %S3 Shear zones
end
OPTVLogs(1).Type = 'All';
OPTVLogs(2).Type = 'Fracture';
OPTVLogs(3).Type = 'S1 Shear-zone';
OPTVLogs(4).Type = 'S3 Shear-zone';

end

