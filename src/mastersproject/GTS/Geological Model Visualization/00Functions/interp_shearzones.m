function [S11_int,S12_int,S13_int,S31_int,S32_int] = interp_shearzones(iv)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% figure(500)
% hold on
figure(iv)
S11_int = dlmread('S11_interp_grid.txt');
S11_plot = surf(S11_int(1:20,:), S11_int(21:40,:),S11_int(41:60,:),'FaceAlpha',0.5,'EdgeAlpha',0.5);
set(S11_plot,'FaceColor','r')
S12_int = dlmread('S12_interp_grid.txt');
S12_plot = surf(S12_int(1:20,:), S12_int(21:40,:),S12_int(41:60,:),'FaceAlpha',0.5,'EdgeAlpha',0.5);
set(S12_plot,'FaceColor','r')
S13_int = dlmread('S13_interp_grid.txt');
S13_plot = surf(S13_int(1:20,:), S13_int(21:40,:),S13_int(41:60,:),'FaceAlpha',0.5,'EdgeAlpha',0.5);
set(S13_plot,'FaceColor','r')
S31_int = dlmread('S31_interp_grid.txt');
S31_plot = surf(S31_int(1:20,:), S31_int(21:40,:),S31_int(41:60,:),'FaceAlpha',0.5,'EdgeAlpha',0.5); 
set(S31_plot,'FaceColor','g')
S32_int = dlmread('S32_interp_grid.txt');
S32_plot = surf(S32_int(1:20,:), S32_int(21:40,:),S32_int(41:60,:),'FaceAlpha',0.5,'EdgeAlpha',0.5);
set(S32_plot,'FaceColor','g')
zlim([0 50])
end

