function [ XX, YY, ZZ, vel ] = Plot2D_seismics(iv, GTS_coordinates)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% warning off

% addpath('..\07_Constraints\01-Geophysics\Seismic_tomography\seismic_2d_aniso\')

xshift = 667410;
yshift = 158880;
zshift = 1730;


eval('aniso_indiv');
load inv2dm_all
% close all
load cov
%load model_res1.mat

% w = ones(size(R,1),1);
% w(d.grid.mx*d.grid.mz+1:end) = 1/4700;
% ww = w*(1./w');
% 
% RR = ww.*R;

%res = full(diag(R));
n_mod = 84; % In this case n_mod is the shift in north-direction
e_mod = 10;
z_mod = 34;
for it=10
    
    vel = reshape(1./load(sprintf('aniso_v%i.txt',it)),d.grid.mx,d.grid.mz);
%     vel = vel/1000;
    vel = vel';
    figure(iv); %set(gca,'fontsize',18);
    [XX, YY]=meshgrid(d.grid.dx+e_mod:d.grid.dx:d.grid.dx*d.grid.mx+e_mod,d.grid.dz+n_mod:d.grid.dz:d.grid.dz*d.grid.mz+n_mod);
    ZZ = ones(size(XX))*z_mod;
    surf(XX,YY,ZZ,vel);
     shading interp;
    caxis([4500 4800])
    
    
    
    tunnel( GTS_coordinates,iv);
    interp_shearzones(iv);
end

