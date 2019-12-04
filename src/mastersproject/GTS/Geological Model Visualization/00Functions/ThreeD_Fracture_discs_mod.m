function [OPTVLogs] = ThreeD_Fracture_discs_mod( OPTVLogs, FBS, SBH, INJ, GEO, PRP)
% function [ output_args ] = ThreeD_Fracture_discs( OPTVLogs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% This function uses the OPTV logs and plots discs into figure(iv)

% Defines figure number where to plot in
for iis = 2:4
%% FBS ====================================================================
OPTVLogs(iis).FBS1(:,4:6) = [1.*sind(OPTVLogs(iis).FBS1(:,2)).*sind(OPTVLogs(iis).FBS1(:,3)), 1.*cosd(OPTVLogs(iis).FBS1(:,2)).*sind(OPTVLogs(iis).FBS1(:,3)),    1.*cosd(OPTVLogs(iis).FBS1(:,3))];
OPTVLogs(iis).FBS2(:,4:6) = [1.*sind(OPTVLogs(iis).FBS2(:,2)).*sind(OPTVLogs(iis).FBS2(:,3)), 1.*cosd(OPTVLogs(iis).FBS2(:,2)).*sind(OPTVLogs(iis).FBS2(:,3)),    1.*cosd(OPTVLogs(iis).FBS2(:,3))];
OPTVLogs(iis).FBS3(:,4:6) = [1.*sind(OPTVLogs(iis).FBS3(:,2)).*sind(OPTVLogs(iis).FBS3(:,3)), 1.*cosd(OPTVLogs(iis).FBS3(:,2)).*sind(OPTVLogs(iis).FBS3(:,3)),    1.*cosd(OPTVLogs(iis).FBS3(:,3))];
OPTVLogs(iis).FBS1(:,7:9) = [FBS(1,1)+OPTVLogs(iis).FBS1(:,1).*sind(FBS(1,6)).*cosd(FBS(1,7)), FBS(1,2)+OPTVLogs(iis).FBS1(:,1).*cosd(FBS(1,6))*cosd(FBS(1,7)), FBS(1,3)+OPTVLogs(iis).FBS1(:,1).*sind(FBS(1,7))];
OPTVLogs(iis).FBS2(:,7:9) = [FBS(2,1)+OPTVLogs(iis).FBS2(:,1).*sind(FBS(2,6)).*cosd(FBS(2,7)), FBS(2,2)+OPTVLogs(iis).FBS2(:,1).*cosd(FBS(2,6))*cosd(FBS(2,7)), FBS(2,3)+OPTVLogs(iis).FBS2(:,1).*sind(FBS(2,7))];
OPTVLogs(iis).FBS3(:,7:9) = [FBS(3,1)+OPTVLogs(iis).FBS3(:,1).*sind(FBS(3,6)).*cosd(FBS(3,7)), FBS(3,2)+OPTVLogs(iis).FBS3(:,1).*cosd(FBS(3,6))*cosd(FBS(3,7)), FBS(3,3)+OPTVLogs(iis).FBS3(:,1).*sind(FBS(3,7))];

%% SBH ====================================================================
OPTVLogs(iis).SBH1(:,4:6) = [1.*sind(OPTVLogs(iis).SBH1(:,2)).*sind(OPTVLogs(iis).SBH1(:,3)), 1.*cosd(OPTVLogs(iis).SBH1(:,2)).*sind(OPTVLogs(iis).SBH1(:,3)),    1.*cosd(OPTVLogs(iis).SBH1(:,3))];
OPTVLogs(iis).SBH3(:,4:6) = [1.*sind(OPTVLogs(iis).SBH3(:,2)).*sind(OPTVLogs(iis).SBH3(:,3)), 1.*cosd(OPTVLogs(iis).SBH3(:,2)).*sind(OPTVLogs(iis).SBH3(:,3)),    1.*cosd(OPTVLogs(iis).SBH3(:,3))];
OPTVLogs(iis).SBH4(:,4:6) = [1.*sind(OPTVLogs(iis).SBH4(:,2)).*sind(OPTVLogs(iis).SBH4(:,3)), 1.*cosd(OPTVLogs(iis).SBH4(:,2)).*sind(OPTVLogs(iis).SBH4(:,3)),    1.*cosd(OPTVLogs(iis).SBH4(:,3))];
OPTVLogs(iis).SBH1(:,7:9) = [SBH(1,1)+OPTVLogs(iis).SBH1(:,1).*sind(SBH(1,6)).*cosd(SBH(1,7)), SBH(1,2)+OPTVLogs(iis).SBH1(:,1).*cosd(SBH(1,6))*cosd(SBH(1,7)), SBH(1,3)+OPTVLogs(iis).SBH1(:,1).*sind(SBH(1,7))];
OPTVLogs(iis).SBH3(:,7:9) = [SBH(2,1)+OPTVLogs(iis).SBH3(:,1).*sind(SBH(2,6)).*cosd(SBH(2,7)), SBH(2,2)+OPTVLogs(iis).SBH3(:,1).*cosd(SBH(2,6))*cosd(SBH(2,7)), SBH(2,3)+OPTVLogs(iis).SBH3(:,1).*sind(SBH(2,7))];
OPTVLogs(iis).SBH4(:,7:9) = [SBH(3,1)+OPTVLogs(iis).SBH4(:,1).*sind(SBH(3,6)).*cosd(SBH(3,7)), SBH(3,2)+OPTVLogs(iis).SBH4(:,1).*cosd(SBH(3,6))*cosd(SBH(3,7)), SBH(3,3)+OPTVLogs(iis).SBH4(:,1).*sind(SBH(3,7))];

%% PRP ====================================================================
OPTVLogs(iis).PRP1(:,4:6) = [1.*sind(OPTVLogs(iis).PRP1(:,2)).*sind(OPTVLogs(iis).PRP1(:,3)), 1.*cosd(OPTVLogs(iis).PRP1(:,2)).*sind(OPTVLogs(iis).PRP1(:,3)),    1.*cosd(OPTVLogs(iis).PRP1(:,3))];
OPTVLogs(iis).PRP2(:,4:6) = [1.*sind(OPTVLogs(iis).PRP2(:,2)).*sind(OPTVLogs(iis).PRP2(:,3)), 1.*cosd(OPTVLogs(iis).PRP2(:,2)).*sind(OPTVLogs(iis).PRP2(:,3)),    1.*cosd(OPTVLogs(iis).PRP2(:,3))];
OPTVLogs(iis).PRP3(:,4:6) = [1.*sind(OPTVLogs(iis).PRP3(:,2)).*sind(OPTVLogs(iis).PRP3(:,3)), 1.*cosd(OPTVLogs(iis).PRP3(:,2)).*sind(OPTVLogs(iis).PRP3(:,3)),    1.*cosd(OPTVLogs(iis).PRP3(:,3))];
OPTVLogs(iis).PRP1(:,7:9) = [PRP(1,1)+OPTVLogs(iis).PRP1(:,1).*sind(PRP(1,6)).*cosd(PRP(1,7)), PRP(1,2)+OPTVLogs(iis).PRP1(:,1).*cosd(PRP(1,6))*cosd(PRP(1,7)), PRP(1,3)+OPTVLogs(iis).PRP1(:,1).*sind(PRP(1,7))];
OPTVLogs(iis).PRP2(:,7:9) = [PRP(2,1)+OPTVLogs(iis).PRP2(:,1).*sind(PRP(2,6)).*cosd(PRP(2,7)), PRP(2,2)+OPTVLogs(iis).PRP2(:,1).*cosd(PRP(2,6))*cosd(PRP(2,7)), PRP(2,3)+OPTVLogs(iis).PRP2(:,1).*sind(PRP(2,7))];
OPTVLogs(iis).PRP3(:,7:9) = [PRP(3,1)+OPTVLogs(iis).PRP3(:,1).*sind(PRP(3,6)).*cosd(PRP(3,7)), PRP(3,2)+OPTVLogs(iis).PRP3(:,1).*cosd(PRP(3,6))*cosd(PRP(3,7)), PRP(3,3)+OPTVLogs(iis).PRP3(:,1).*sind(PRP(3,7))];



%% INJ ====================================================================
OPTVLogs(iis).INJ1(:,4:6) = [1.*sind(OPTVLogs(iis).INJ1(:,2)).*sind(OPTVLogs(iis).INJ1(:,3)), 1.*cosd(OPTVLogs(iis).INJ1(:,2)).*sind(OPTVLogs(iis).INJ1(:,3)),    1.*cosd(OPTVLogs(iis).INJ1(:,3))];
OPTVLogs(iis).INJ2(:,4:6) = [1.*sind(OPTVLogs(iis).INJ2(:,2)).*sind(OPTVLogs(iis).INJ2(:,3)), 1.*cosd(OPTVLogs(iis).INJ2(:,2)).*sind(OPTVLogs(iis).INJ2(:,3)),    1.*cosd(OPTVLogs(iis).INJ2(:,3))];

OPTVLogs(iis).INJ1(:,7:9) = [INJ(1,1)+OPTVLogs(iis).INJ1(:,1).*sind(INJ(1,6)).*cosd(INJ(1,7)), INJ(1,2)+OPTVLogs(iis).INJ1(:,1).*cosd(INJ(1,6))*cosd(INJ(1,7)), INJ(1,3)+OPTVLogs(iis).INJ1(:,1).*sind(INJ(1,7))];
OPTVLogs(iis).INJ2(:,7:9) = [INJ(2,1)+OPTVLogs(iis).INJ2(:,1).*sind(INJ(2,6)).*cosd(INJ(2,7)), INJ(2,2)+OPTVLogs(iis).INJ2(:,1).*cosd(INJ(2,6))*cosd(INJ(2,7)), INJ(2,3)+OPTVLogs(iis).INJ2(:,1).*sind(INJ(2,7))];


%% GEO ====================================================================
OPTVLogs(iis).GEO1(:,4:6) = [1.*sind(OPTVLogs(iis).GEO1(:,2)).*sind(OPTVLogs(iis).GEO1(:,3)), 1.*cosd(OPTVLogs(iis).GEO1(:,2)).*sind(OPTVLogs(iis).GEO1(:,3)),    1.*cosd(OPTVLogs(iis).GEO1(:,3))];
OPTVLogs(iis).GEO2(:,4:6) = [1.*sind(OPTVLogs(iis).GEO2(:,2)).*sind(OPTVLogs(iis).GEO2(:,3)), 1.*cosd(OPTVLogs(iis).GEO2(:,2)).*sind(OPTVLogs(iis).GEO2(:,3)),    1.*cosd(OPTVLogs(iis).GEO2(:,3))];
OPTVLogs(iis).GEO3(:,4:6) = [1.*sind(OPTVLogs(iis).GEO3(:,2)).*sind(OPTVLogs(iis).GEO3(:,3)), 1.*cosd(OPTVLogs(iis).GEO3(:,2)).*sind(OPTVLogs(iis).GEO3(:,3)),    1.*cosd(OPTVLogs(iis).GEO3(:,3))];
OPTVLogs(iis).GEO4(:,4:6) = [1.*sind(OPTVLogs(iis).GEO4(:,2)).*sind(OPTVLogs(iis).GEO4(:,3)), 1.*cosd(OPTVLogs(iis).GEO4(:,2)).*sind(OPTVLogs(iis).GEO4(:,3)),    1.*cosd(OPTVLogs(iis).GEO4(:,3))];
OPTVLogs(iis).GEO1(:,7:9) = [GEO(1,1)+OPTVLogs(iis).GEO1(:,1).*sind(GEO(1,6)).*cosd(GEO(1,7)), GEO(1,2)+OPTVLogs(iis).GEO1(:,1).*cosd(GEO(1,6))*cosd(GEO(1,7)), GEO(1,3)+OPTVLogs(iis).GEO1(:,1).*sind(GEO(1,7))];
OPTVLogs(iis).GEO2(:,7:9) = [GEO(2,1)+OPTVLogs(iis).GEO2(:,1).*sind(GEO(2,6)).*cosd(GEO(2,7)), GEO(2,2)+OPTVLogs(iis).GEO2(:,1).*cosd(GEO(2,6))*cosd(GEO(2,7)), GEO(2,3)+OPTVLogs(iis).GEO2(:,1).*sind(GEO(2,7))];
OPTVLogs(iis).GEO3(:,7:9) = [GEO(3,1)+OPTVLogs(iis).GEO3(:,1).*sind(GEO(3,6)).*cosd(GEO(3,7)), GEO(3,2)+OPTVLogs(iis).GEO3(:,1).*cosd(GEO(3,6))*cosd(GEO(3,7)), GEO(3,3)+OPTVLogs(iis).GEO3(:,1).*sind(GEO(3,7))];
OPTVLogs(iis).GEO4(:,7:9) = [GEO(4,1)+OPTVLogs(iis).GEO4(:,1).*sind(GEO(4,6)).*cosd(GEO(4,7)), GEO(4,2)+OPTVLogs(iis).GEO4(:,1).*cosd(GEO(4,6))*cosd(GEO(4,7)), GEO(4,3)+OPTVLogs(iis).GEO4(:,1).*sind(GEO(4,7))];


end
end

