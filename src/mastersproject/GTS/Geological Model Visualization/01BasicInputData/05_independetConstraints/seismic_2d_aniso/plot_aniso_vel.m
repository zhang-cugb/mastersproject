% plot seismic velocity tomogram

clear all

xshift = 667410;
yshift = 158880;
zshift = 1730;


eval('aniso_indiv');
load inv2dm_all
close all
load cov
%load model_res1.mat

% w = ones(size(R,1),1);
% w(d.grid.mx*d.grid.mz+1:end) = 1/4700;
% ww = w*(1./w');
% 
% RR = ww.*R;

%res = full(diag(R));

for it=10
    
    vel = reshape(1./load(sprintf('aniso_v%i.txt',it)),d.grid.mx,d.grid.mz);
    
    figure(1); clf; set(gca,'fontsize',18);
    subplot(2,2,1); hold on;
    imagesc(d.grid.dx:d.grid.dx:d.grid.dx*d.grid.mx,d.grid.dz:d.grid.dz:d.grid.dz*d.grid.mz,vel');
    plot(all.geo.x,all.geo.z,'k.')
    axis equal tight
    patch(x-.5,y,'w','edgecolor','none','facealpha',.7)
    caxis([4500 4800])
    h = colorbar;
    title(sprintf('Velocity, it %i',it))
    xlabel('Easting [m]')
    ylabel('Northing [m]')
    set(gca,'xlim',[1 63],'ylim',[1 62])
    box on
    
    subplot(2,2,2);  hold on;
    eps = reshape(load(sprintf('aniso_e%i.txt',it)),d.grid.mx,d.grid.mz);
    imagesc(d.grid.dx:d.grid.dx:d.grid.dx*d.grid.mx,d.grid.dz:d.grid.dz:d.grid.dz*d.grid.mz,eps');
    plot(all.geo.x,all.geo.z,'k.')
    axis equal tight
    patch(x-.5,y,'w','edgecolor','none','facealpha',.7)
    caxis([0.05 .09])
    h = colorbar;
    title(sprintf('epsilon'))
    xlabel('Easting [m]')
    ylabel('Northing [m]')
    set(gca,'xlim',[1 63],'ylim',[1 62])
    box on
    
    subplot(2,2,3);  hold on;
    delta = reshape(load(sprintf('aniso_d%i.txt',it)),d.grid.mx,d.grid.mz);
    imagesc(d.grid.dx:d.grid.dx:d.grid.dx*d.grid.mx,d.grid.dz:d.grid.dz:d.grid.dz*d.grid.mz,delta');
    plot(all.geo.x,all.geo.z,'k.')
    axis equal tight
    patch(x-.5,y,'w','edgecolor','none','facealpha',.7)
    %    caxis([0.07 .1])
    h = colorbar;
    title(sprintf('delta'))
    xlabel('Easting [m]')
    ylabel('Northing [m]')
    set(gca,'xlim',[1 63],'ylim',[1 62])
    box on
    
    subplot(2,2,4)
    st = load(sprintf('aniso_st%i.txt',it));
    plot(st*1000,'ko');
    title(sprintf('Source statics'))
    print -dpng -r300 9_dir_st_corr
    %    pause()
    %         colormap jet
    %          subplot(2,2,1);
    %          caxis([4400 5000])
    %          subplot(2,2,2);
    %          caxis([.04 .085])
    %          subplot(2,2,3);
    %          caxis([.08 .14])
    %          subplot(2,2,4);
    %          caxis([120 145])
    %         print -dpng -r300 15_joccam_all
    
    
end

