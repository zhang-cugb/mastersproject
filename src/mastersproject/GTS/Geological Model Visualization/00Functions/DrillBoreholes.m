function [ FBS, INJ, GEO, PRP, SBH] = DrillBoreholes( BH_import, BH_plot,GTS_coordinates,num_fig )
% This function imports borehole information and plots the chosen boreholes

% Import Borehole Data
for ii = 1:length(BH_import)
load(strcat(BH_import{ii},'.txt'));
end

% Plot Boreholes
for ii = 1:length(BH_plot)
    a = load(strcat(BH_plot{ii},'.txt'));
a(:,1) = a(:,1) - GTS_coordinates.x;
a(:,2) = a(:,2) - GTS_coordinates.y;
a(:,3) = a(:,3) - GTS_coordinates.z;

% Define Colorcode
if strcmp(BH_plot(ii),{'FBS'})
colorcode = 'r';
elseif strcmp(BH_plot(ii),{'SBH'})
    colorcode = 'k';
elseif strcmp(BH_plot(ii),{'INJ'})
%     colorcode = [0.3 0.3 0.3];
colorcode = 'g';
elseif strcmp(BH_plot(ii),{'PRP'})
    colorcode = 'm';
elseif strcmp(BH_plot(ii),{'GEO'})
    colorcode = 'y';
end

    for iii = 1:num_fig
       figure(iii)
       for iv = 1:length(a(:,1))
              plot3([a(iv,1) a(iv,1)+a(iv,4).*sind(a(iv,6)).*cosd(a(iv,7))],[a(iv,2)  a(iv,2)+a(iv,4).*cosd(a(iv,6))*cosd(a(iv,7))],[a(iv,3) a(iv,3)+a(iv,4)*sind(a(iv,7))],'Color',colorcode,'Linewidth',2)
       end
    end
end


FBS(:,1) = FBS(:,1) - GTS_coordinates.x;
FBS(:,2) = FBS(:,2) - GTS_coordinates.y;
FBS(:,3) = FBS(:,3) - GTS_coordinates.z;

GEO(:,1) = GEO(:,1) - GTS_coordinates.x;
GEO(:,2) = GEO(:,2) - GTS_coordinates.y;
GEO(:,3) = GEO(:,3) - GTS_coordinates.z;

PRP(:,1) = PRP(:,1) - GTS_coordinates.x;
PRP(:,2) = PRP(:,2) - GTS_coordinates.y;
PRP(:,3) = PRP(:,3) - GTS_coordinates.z;

INJ(:,1) = INJ(:,1) - GTS_coordinates.x;
INJ(:,2) = INJ(:,2) - GTS_coordinates.y;
INJ(:,3) = INJ(:,3) - GTS_coordinates.z;

SBH(:,1) = SBH(:,1) - GTS_coordinates.x;
SBH(:,2) = SBH(:,2) - GTS_coordinates.y;
SBH(:,3) = SBH(:,3) - GTS_coordinates.z;


end

