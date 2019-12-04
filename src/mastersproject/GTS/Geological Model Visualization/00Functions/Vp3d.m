function [ FracturePerMeter ] = frac_density3d( FBS, SBH, INJ, GEO, PRP ,iv )

% Import Fracture per Meter log
[ FracturePerMeter ] = importFracDensity;

figure(iv)
addpath('..\07_Constraints\01-Geophysics\Seismic_tomography\3D_seismic')
addpath('..\07_Constraints\01-Geophysics\Seismic_tomography\3D_seismic\seismic')

% plot 3D seismic
[Vp ,x_seis,y_seis,z_seis] = occ_show('tomo3d',10,68,108,[5]);
caxis([4.6 5.4])



