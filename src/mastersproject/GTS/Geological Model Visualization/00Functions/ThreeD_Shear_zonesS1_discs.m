function [OPTVLogs] = ThreeD_Shear_zonesS1_discs( OPTVLogs, FBS, SBH, INJ, GEO, PRP,iv )
% function [ output_args ] = ThreeD_Fracture_discs( OPTVLogs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% This function uses the OPTV logs and plots discs into 
radius= 2;
figure(iv)
%% FBS ====================================================================
OPTVLogs(3).FBS1(:,4:6) = [1.*sind(OPTVLogs(3).FBS1(:,2)).*sind(OPTVLogs(3).FBS1(:,3)), 1.*cosd(OPTVLogs(3).FBS1(:,2)).*sind(OPTVLogs(3).FBS1(:,3)),    1.*cosd(OPTVLogs(3).FBS1(:,3))];
OPTVLogs(3).FBS2(:,4:6) = [1.*sind(OPTVLogs(3).FBS2(:,2)).*sind(OPTVLogs(3).FBS2(:,3)), 1.*cosd(OPTVLogs(3).FBS2(:,2)).*sind(OPTVLogs(3).FBS2(:,3)),    1.*cosd(OPTVLogs(3).FBS2(:,3))];
OPTVLogs(3).FBS3(:,4:6) = [1.*sind(OPTVLogs(3).FBS3(:,2)).*sind(OPTVLogs(3).FBS3(:,3)), 1.*cosd(OPTVLogs(3).FBS3(:,2)).*sind(OPTVLogs(3).FBS3(:,3)),    1.*cosd(OPTVLogs(3).FBS3(:,3))];
OPTVLogs(3).FBS1(:,7:9) = [FBS(1,1)+OPTVLogs(3).FBS1(:,1).*sind(FBS(1,6)).*cosd(FBS(1,7)), FBS(1,2)+OPTVLogs(3).FBS1(:,1).*cosd(FBS(1,6))*cosd(FBS(1,7)), FBS(1,3)+OPTVLogs(3).FBS1(:,1).*sind(FBS(1,7))];
OPTVLogs(3).FBS2(:,7:9) = [FBS(2,1)+OPTVLogs(3).FBS2(:,1).*sind(FBS(2,6)).*cosd(FBS(2,7)), FBS(2,2)+OPTVLogs(3).FBS2(:,1).*cosd(FBS(2,6))*cosd(FBS(2,7)), FBS(2,3)+OPTVLogs(3).FBS2(:,1).*sind(FBS(2,7))];
OPTVLogs(3).FBS3(:,7:9) = [FBS(3,1)+OPTVLogs(3).FBS3(:,1).*sind(FBS(3,6)).*cosd(FBS(3,7)), FBS(3,2)+OPTVLogs(3).FBS3(:,1).*cosd(FBS(3,6))*cosd(FBS(3,7)), FBS(3,3)+OPTVLogs(3).FBS3(:,1).*sind(FBS(3,7))];
for ii = 1:length(OPTVLogs(3).FBS1(:,7))

center= [OPTVLogs(3).FBS1(ii,7) OPTVLogs(3).FBS1(ii,8), OPTVLogs(3).FBS1(ii,9)];
normal = [OPTVLogs(3).FBS1(ii,4), OPTVLogs(3).FBS1(ii,5), OPTVLogs(3).FBS1(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


for ii = 1:length(OPTVLogs(3).FBS2(:,7))

center= [OPTVLogs(3).FBS2(ii,7) OPTVLogs(3).FBS2(ii,8), OPTVLogs(3).FBS2(ii,9)];
normal = [OPTVLogs(3).FBS1(ii,4), OPTVLogs(3).FBS2(ii,5), OPTVLogs(3).FBS2(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end

for ii = 1:length(OPTVLogs(3).FBS3(:,7))

center= [OPTVLogs(3).FBS3(ii,7) OPTVLogs(3).FBS3(ii,8), OPTVLogs(3).FBS3(ii,9)];
normal = [OPTVLogs(3).FBS3(ii,4), OPTVLogs(3).FBS3(ii,5), OPTVLogs(3).FBS3(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


%% SBH ====================================================================
OPTVLogs(3).SBH1(:,4:6) = [1.*sind(OPTVLogs(3).SBH1(:,2)).*sind(OPTVLogs(3).SBH1(:,3)), 1.*cosd(OPTVLogs(3).SBH1(:,2)).*sind(OPTVLogs(3).SBH1(:,3)),    1.*cosd(OPTVLogs(3).SBH1(:,3))];
OPTVLogs(3).SBH3(:,4:6) = [1.*sind(OPTVLogs(3).SBH3(:,2)).*sind(OPTVLogs(3).SBH3(:,3)), 1.*cosd(OPTVLogs(3).SBH3(:,2)).*sind(OPTVLogs(3).SBH3(:,3)),    1.*cosd(OPTVLogs(3).SBH3(:,3))];
OPTVLogs(3).SBH4(:,4:6) = [1.*sind(OPTVLogs(3).SBH4(:,2)).*sind(OPTVLogs(3).SBH4(:,3)), 1.*cosd(OPTVLogs(3).SBH4(:,2)).*sind(OPTVLogs(3).SBH4(:,3)),    1.*cosd(OPTVLogs(3).SBH4(:,3))];
OPTVLogs(3).SBH1(:,7:9) = [SBH(1,1)+OPTVLogs(3).SBH1(:,1).*sind(SBH(1,6)).*cosd(SBH(1,7)), SBH(1,2)+OPTVLogs(3).SBH1(:,1).*cosd(SBH(1,6))*cosd(SBH(1,7)), SBH(1,3)+OPTVLogs(3).SBH1(:,1).*sind(SBH(1,7))];
OPTVLogs(3).SBH3(:,7:9) = [SBH(2,1)+OPTVLogs(3).SBH3(:,1).*sind(SBH(2,6)).*cosd(SBH(2,7)), SBH(2,2)+OPTVLogs(3).SBH3(:,1).*cosd(SBH(2,6))*cosd(SBH(2,7)), SBH(2,3)+OPTVLogs(3).SBH3(:,1).*sind(SBH(2,7))];
OPTVLogs(3).SBH4(:,7:9) = [SBH(3,1)+OPTVLogs(3).SBH4(:,1).*sind(SBH(3,6)).*cosd(SBH(3,7)), SBH(3,2)+OPTVLogs(3).SBH4(:,1).*cosd(SBH(3,6))*cosd(SBH(3,7)), SBH(3,3)+OPTVLogs(3).SBH4(:,1).*sind(SBH(3,7))];
for ii = 1:length(OPTVLogs(3).SBH1(:,7))

center= [OPTVLogs(3).SBH1(ii,7) OPTVLogs(3).SBH1(ii,8), OPTVLogs(3).SBH1(ii,9)];
normal = [OPTVLogs(3).SBH1(ii,4), OPTVLogs(3).SBH1(ii,5), OPTVLogs(3).SBH1(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


for ii = 1:length(OPTVLogs(3).SBH4(:,7))

center= [OPTVLogs(3).SBH4(ii,7) OPTVLogs(3).SBH4(ii,8), OPTVLogs(3).SBH4(ii,9)];
normal = [OPTVLogs(3).SBH4(ii,4), OPTVLogs(3).SBH4(ii,5), OPTVLogs(3).SBH4(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end

for ii = 1:length(OPTVLogs(3).SBH4(:,7))

center= [OPTVLogs(3).SBH4(ii,7) OPTVLogs(3).SBH4(ii,8), OPTVLogs(3).SBH4(ii,9)];
normal = [OPTVLogs(3).SBH4(ii,4), OPTVLogs(3).SBH4(ii,5), OPTVLogs(3).SBH4(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end

%% PRP ====================================================================
OPTVLogs(3).PRP1(:,4:6) = [1.*sind(OPTVLogs(3).PRP1(:,2)).*sind(OPTVLogs(3).PRP1(:,3)), 1.*cosd(OPTVLogs(3).PRP1(:,2)).*sind(OPTVLogs(3).PRP1(:,3)),    1.*cosd(OPTVLogs(3).PRP1(:,3))];
OPTVLogs(3).PRP2(:,4:6) = [1.*sind(OPTVLogs(3).PRP2(:,2)).*sind(OPTVLogs(3).PRP2(:,3)), 1.*cosd(OPTVLogs(3).PRP2(:,2)).*sind(OPTVLogs(3).PRP2(:,3)),    1.*cosd(OPTVLogs(3).PRP2(:,3))];
OPTVLogs(3).PRP3(:,4:6) = [1.*sind(OPTVLogs(3).PRP3(:,2)).*sind(OPTVLogs(3).PRP3(:,3)), 1.*cosd(OPTVLogs(3).PRP3(:,2)).*sind(OPTVLogs(3).PRP3(:,3)),    1.*cosd(OPTVLogs(3).PRP3(:,3))];
OPTVLogs(3).PRP1(:,7:9) = [PRP(1,1)+OPTVLogs(3).PRP1(:,1).*sind(PRP(1,6)).*cosd(PRP(1,7)), PRP(1,2)+OPTVLogs(3).PRP1(:,1).*cosd(PRP(1,6))*cosd(PRP(1,7)), PRP(1,3)+OPTVLogs(3).PRP1(:,1).*sind(PRP(1,7))];
OPTVLogs(3).PRP2(:,7:9) = [PRP(2,1)+OPTVLogs(3).PRP2(:,1).*sind(PRP(2,6)).*cosd(PRP(2,7)), PRP(2,2)+OPTVLogs(3).PRP2(:,1).*cosd(PRP(2,6))*cosd(PRP(2,7)), PRP(2,3)+OPTVLogs(3).PRP2(:,1).*sind(PRP(2,7))];
OPTVLogs(3).PRP3(:,7:9) = [PRP(3,1)+OPTVLogs(3).PRP3(:,1).*sind(PRP(3,6)).*cosd(PRP(3,7)), PRP(3,2)+OPTVLogs(3).PRP3(:,1).*cosd(PRP(3,6))*cosd(PRP(3,7)), PRP(3,3)+OPTVLogs(3).PRP3(:,1).*sind(PRP(3,7))];
for ii = 1:length(OPTVLogs(3).PRP1(:,7))

center= [OPTVLogs(3).PRP1(ii,7) OPTVLogs(3).PRP1(ii,8), OPTVLogs(3).PRP1(ii,9)];
normal = [OPTVLogs(3).PRP1(ii,4), OPTVLogs(3).PRP1(ii,5), OPTVLogs(3).PRP1(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


for ii = 1:length(OPTVLogs(3).PRP2(:,7))

center= [OPTVLogs(3).PRP2(ii,7) OPTVLogs(3).PRP2(ii,8), OPTVLogs(3).PRP2(ii,9)];
normal = [OPTVLogs(3).PRP2(ii,4), OPTVLogs(3).PRP2(ii,5), OPTVLogs(3).PRP2(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end

for ii = 1:length(OPTVLogs(3).PRP3(:,7))

center= [OPTVLogs(3).PRP3(ii,7) OPTVLogs(3).PRP3(ii,8), OPTVLogs(3).PRP3(ii,9)];
normal = [OPTVLogs(3).PRP3(ii,4), OPTVLogs(3).PRP3(ii,5), OPTVLogs(3).PRP3(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


%% INJ ====================================================================
OPTVLogs(3).INJ1(:,4:6) = [1.*sind(OPTVLogs(3).INJ1(:,2)).*sind(OPTVLogs(3).INJ1(:,3)), 1.*cosd(OPTVLogs(3).INJ1(:,2)).*sind(OPTVLogs(3).INJ1(:,3)),    1.*cosd(OPTVLogs(3).INJ1(:,3))];
OPTVLogs(3).INJ2(:,4:6) = [1.*sind(OPTVLogs(3).INJ2(:,2)).*sind(OPTVLogs(3).INJ2(:,3)), 1.*cosd(OPTVLogs(3).INJ2(:,2)).*sind(OPTVLogs(3).INJ2(:,3)),    1.*cosd(OPTVLogs(3).INJ2(:,3))];

OPTVLogs(3).INJ1(:,7:9) = [INJ(1,1)+OPTVLogs(3).INJ1(:,1).*sind(INJ(1,6)).*cosd(INJ(1,7)), INJ(1,2)+OPTVLogs(3).INJ1(:,1).*cosd(INJ(1,6))*cosd(INJ(1,7)), INJ(1,3)+OPTVLogs(3).INJ1(:,1).*sind(INJ(1,7))];
OPTVLogs(3).INJ2(:,7:9) = [INJ(2,1)+OPTVLogs(3).INJ2(:,1).*sind(INJ(2,6)).*cosd(INJ(2,7)), INJ(2,2)+OPTVLogs(3).INJ2(:,1).*cosd(INJ(2,6))*cosd(INJ(2,7)), INJ(2,3)+OPTVLogs(3).INJ2(:,1).*sind(INJ(2,7))];

for ii = 1:length(OPTVLogs(3).INJ1(:,7))

center= [OPTVLogs(3).INJ1(ii,7) OPTVLogs(3).INJ1(ii,8), OPTVLogs(3).INJ1(ii,9)];
normal = [OPTVLogs(3).INJ1(ii,4), OPTVLogs(3).INJ1(ii,5), OPTVLogs(3).INJ1(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


for ii = 1:length(OPTVLogs(3).INJ2(:,7))

center= [OPTVLogs(3).INJ2(ii,7) OPTVLogs(3).INJ2(ii,8), OPTVLogs(3).INJ2(ii,9)];
normal = [OPTVLogs(3).INJ2(ii,4), OPTVLogs(3).INJ2(ii,5), OPTVLogs(3).INJ2(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


%% GEO ====================================================================
OPTVLogs(3).GEO1(:,4:6) = [1.*sind(OPTVLogs(3).GEO1(:,2)).*sind(OPTVLogs(3).GEO1(:,3)), 1.*cosd(OPTVLogs(3).GEO1(:,2)).*sind(OPTVLogs(3).GEO1(:,3)),    1.*cosd(OPTVLogs(3).GEO1(:,3))];
OPTVLogs(3).GEO2(:,4:6) = [1.*sind(OPTVLogs(3).GEO2(:,2)).*sind(OPTVLogs(3).GEO2(:,3)), 1.*cosd(OPTVLogs(3).GEO2(:,2)).*sind(OPTVLogs(3).GEO2(:,3)),    1.*cosd(OPTVLogs(3).GEO2(:,3))];
OPTVLogs(3).GEO3(:,4:6) = [1.*sind(OPTVLogs(3).GEO3(:,2)).*sind(OPTVLogs(3).GEO3(:,3)), 1.*cosd(OPTVLogs(3).GEO3(:,2)).*sind(OPTVLogs(3).GEO3(:,3)),    1.*cosd(OPTVLogs(3).GEO3(:,3))];
OPTVLogs(3).GEO4(:,4:6) = [1.*sind(OPTVLogs(3).GEO4(:,2)).*sind(OPTVLogs(3).GEO4(:,3)), 1.*cosd(OPTVLogs(3).GEO4(:,2)).*sind(OPTVLogs(3).GEO4(:,3)),    1.*cosd(OPTVLogs(3).GEO4(:,3))];
OPTVLogs(3).GEO1(:,7:9) = [GEO(1,1)+OPTVLogs(3).GEO1(:,1).*sind(GEO(1,6)).*cosd(GEO(1,7)), GEO(1,2)+OPTVLogs(3).GEO1(:,1).*cosd(GEO(1,6))*cosd(GEO(1,7)), GEO(1,3)+OPTVLogs(3).GEO1(:,1).*sind(GEO(1,7))];
OPTVLogs(3).GEO2(:,7:9) = [GEO(2,1)+OPTVLogs(3).GEO2(:,1).*sind(GEO(2,6)).*cosd(GEO(2,7)), GEO(2,2)+OPTVLogs(3).GEO2(:,1).*cosd(GEO(2,6))*cosd(GEO(2,7)), GEO(2,3)+OPTVLogs(3).GEO2(:,1).*sind(GEO(2,7))];
OPTVLogs(3).GEO3(:,7:9) = [GEO(3,1)+OPTVLogs(3).GEO3(:,1).*sind(GEO(3,6)).*cosd(GEO(3,7)), GEO(3,2)+OPTVLogs(3).GEO3(:,1).*cosd(GEO(3,6))*cosd(GEO(3,7)), GEO(3,3)+OPTVLogs(3).GEO3(:,1).*sind(GEO(3,7))];
OPTVLogs(3).GEO4(:,7:9) = [GEO(4,1)+OPTVLogs(3).GEO4(:,1).*sind(GEO(4,6)).*cosd(GEO(4,7)), GEO(4,2)+OPTVLogs(3).GEO4(:,1).*cosd(GEO(4,6))*cosd(GEO(4,7)), GEO(4,3)+OPTVLogs(3).GEO4(:,1).*sind(GEO(4,7))];
for ii = 1:length(OPTVLogs(3).GEO1(:,7))

center= [OPTVLogs(3).GEO1(ii,7) OPTVLogs(3).GEO1(ii,8), OPTVLogs(3).GEO1(ii,9)];
normal = [OPTVLogs(3).GEO1(ii,4), OPTVLogs(3).GEO1(ii,5), OPTVLogs(3).GEO1(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end


for ii = 1:length(OPTVLogs(3).GEO2(:,7))

center= [OPTVLogs(3).GEO2(ii,7) OPTVLogs(3).GEO2(ii,8), OPTVLogs(3).GEO2(ii,9)];
normal = [OPTVLogs(3).GEO2(ii,4), OPTVLogs(3).GEO2(ii,5), OPTVLogs(3).GEO2(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end

for ii = 1:length(OPTVLogs(3).GEO3(:,7))

center= [OPTVLogs(3).GEO3(ii,7) OPTVLogs(3).GEO3(ii,8), OPTVLogs(3).GEO3(ii,9)];
normal = [OPTVLogs(3).GEO3(ii,4), OPTVLogs(3).GEO3(ii,5), OPTVLogs(3).GEO3(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end

for ii = 1:length(OPTVLogs(3).GEO4(:,7))

center= [OPTVLogs(3).GEO4(ii,7) OPTVLogs(3).GEO4(ii,8), OPTVLogs(3).GEO4(ii,9)];
normal = [OPTVLogs(3).GEO4(ii,4), OPTVLogs(3).GEO4(ii,5), OPTVLogs(3).GEO4(ii,6)];
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));

patch(points(1,:),points(2,:),points(3,:),'r-');
end
end

