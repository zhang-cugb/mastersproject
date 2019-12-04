%% -------- Build geological model GTS ------------------------------------
%
%   Version: March 2018 
%   H. Krietsch
%   The geological model is linked to the publication:
%   " Comprehensive geological data for a fractured crystalline rock mass
%   analog for hydraulic stimulation experiments" Krietsch et al. (in
%   review) -> planned for a submission to Scientific Data - Nature
%
%   The geological raw data are published via ETH research collection. The
%   title of the data set is:
%   "Comprehensive geological dataset for a fractured crystalline rock 
%   volume at the Grimsel Test Site" Krietsch et al. 2018 
%   DOI: 10.3929/ethz-b-000243199 
%
%   If you use figures produced with this visualization tool, or geological
%   data presented in the published data set, please use the journal
%   publication for citation.
%   
%   Contact: hannes.krietsch@erdw.ethz.ch
% 
% 
%--------------------------------------------------------------------------
% clc
clear
close all
warning off


%% preparation ------------------------------------------------------------
% Include paths of basic functions and Input data
addpath('00Functions/')
addpath(genpath('01BasicInputData/'))


% number of figures
num_fig =7;    

for ii = 1:num_fig
figure(ii)
hold on
end

% Origin of coordinate system     -> We use here a new coordinate system
% that can be transfered into swiss grid by adding the number below to the
% coordinates in the figures.
GTS_coordinates.x = 667400;
GTS_coordinates.y = 158800;
GTS_coordinates.z = 1700;

%% ================== Import section ======================================
% Here all required data will be imported

%% 1. Step: Define Tunnel =================================================
% Plot tunnel -> The tunnels will be plotted in all figures.
tunnel( GTS_coordinates,num_fig);

%% 2. Step: Drill Borehole ================================================
% All boreholes have to be imported here in order to use their coordinates
% for further calculations. Don't modify BH_import!!. To decide which
% boreholes will be displayed in the Figure use BH_plot. The Boreholes will
% than plotted in every Figure.
BH_import = [{'FBS'},{'SBH'},{'INJ'},{'PRP'},{'GEO'}];

% Chose boreholes that will be displayed
% BH_plot = [{'INJ'}]; % This line would just plot INJ boreholes
BH_plot = BH_import;

% This function drilles the boreholes
% The colorcode for the borehole groups is define in DrillBoreholes in the
% section "Define Colorcode".
[ FBS, INJ, GEO, PRP, SBH] = DrillBoreholes( BH_import, BH_plot ,GTS_coordinates,num_fig);


%% 3. Step: Import structures from optv logs ==============================
% This function imports mapped geological structures from OPTV logs. By
% default the structures from all boreholes will be imported. To modify
% change "Choose Borehoels to import data" section in function.
% The output table "OPTVLogs" contains all data and additionally separated 
% them with respect to Fractures, S1 shear zones and S3 shear zones.
[OPTVLogs ] = importOPTVstructures;


%% 4. Step: Import data from geodetic mapping along tunnel walls ==========
% The output table includes x/y/z coordinates of the mapped structure, as
% well as information about: Dip-direction/dip, in which tunnel the
% measurement was done, and which shear zone set it belongs to.
[ SZ_tunnel] = Tunnel_intersections(GTS_coordinates);


%% ======================= INFO ===========================================
% From this point on, the user can choose, which figure should show which
% data. For this purpose we introduce the variable "iv".
% =========================================================================

%% 5. Step: Calculate locations of geological features ====================
% In this section the true locations of all mapped geological features are
% calculated. The structures are displayed as discs with true local
% dip-direction and dip of the corresponding structure. By default the 
% discs have a radius of 1 m (fractures) and 2 m (Shear zones. This can be 
% changed in each function. The default colorcode is: fractures (b), 
% S1 shear zones (r), S3 shear zones(g).

% Calculate discs around OPTV Fractures
iv = 2; % Define Figure for plotting Fracture discs
ThreeD_Fracture_discs( OPTVLogs, FBS, SBH, INJ, GEO, PRP, iv);


% Calculate discs around OPTV Shear zones
iv = 2; % Define Figure for plotting Fracture discs
ThreeD_Shear_zonesS1_discs( OPTVLogs, FBS, SBH, INJ, GEO, PRP, iv);
ThreeD_Shear_zonesS3_discs( OPTVLogs, FBS, SBH, INJ, GEO, PRP, iv);
 
% Calculate discs around structures from tunnels.
iv = 2; % Define Figure for plotting Fracture discs
[SZ_tunnel] = ThreeD_Shear_zones_tunnels(SZ_tunnel,iv);


%% 6. Step: First linear interpolation between geological observations ====
% These functions define linear patches between mapped shear zone
% coordinates that belong to the same set. The default colorcode is: 
% S1 shear zones (r), S3 shear zones (g). NOTE: this visualization
% represents planes that goes throught the middle of the shear zones.
% The visualization does not contain information about true thickness of
% the shear zones!! Also, it does not take true local orientations into
% account.
iv = 3; % Define Figure for plotting Shear zone patches
% Plot S1 shear zone patches
[S11, S12, S13] =S1_shearzones_patches(SZ_tunnel, FBS, SBH, PRP, INJ, GEO,iv);

% Plot S3 shear zone patches
[S31, S32] =S3_shearzones_patches(SZ_tunnel, FBS, SBH, PRP, INJ, GEO,iv);

%% 7. Step: Plot interpolated shear zones based on new coordinates ========
% This interpolation includes the true locate orientations of the
% structures and uses third order polynomial functions to interpolate
% between the coordinates. This produces new coordinates within the volume
% for all shear zones.

for iv = 4 :7 % Define Figure -> if you want to plot into various figures, use a loop.
[S11_int,S12_int,S13_int,S31_int,S32_int] = interp_shearzones(iv);
end

%% 8. Step: Plot Fracture Densities
% The data included in this density plot are based on core and OPTV logs.
% The fracture densities are plotted as scatter plots
iv = 4;     % Define Figure
[ FracturePerMeter ] = frac_density3d( FBS, SBH, INJ, GEO, PRP ,iv );

%% 9. Step: Plot Seismc 2D tomography into figure
% The presented data were obtained during a tunnel-to-tunnel seismic
% tomography, conducted form AU to VE tunnel.
iv =5;

Plot2D_seismics(iv,GTS_coordinates);
% frac_density( FBS, SBH, INJ, GEO, PRP ,iv );

%% 10. Step: Include 3D Vp tomography
% These data were obtained from borehole sparker data
iv = 6;
figure(iv)
[Vp ,x_seis,y_seis,z_seis] = occ_show3D('tomo3d',10,68,108,5);
caxis([4.6 5.4]);

%% 11. Step: Plot GPR IMAGE obtained in AU tunnel
ivGPR = 7;
plot_GPR(ivGPR);


%% Plot Modification=======================================================
% All plots will be presented the same way, except for GPR image
for ii =1:num_fig
figure(ii)
axis equal
grid on
grid minor
xlabel('Easting +667400')
ylabel('Northing +158800')
zlabel('Height +1700')
% material shiny
light;
axis([00 80 55 150 0 40]);
view(90-25,30)
end

% ------ GPR Image---------------
figure(ivGPR)
title('GPR data and interpolated shear zones')
axis equal
grid on
grid minor
xlabel('Easting +667400')
ylabel('Northing +158800')
zlabel('Height +1700')
% material shiny
light;
axis([50 80 55 150 0 40]);
view(-90,45);
camlight;
