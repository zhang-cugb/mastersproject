function [Vp ,x_seis,y_seis,z_seis] = occ_show3D(pfile,it,xp,yp,zp,iso)
% function plotWiden(modelvec,data_type)
% modelvec should contain the data to be plotted
% data_type is a string with on of
% data_type: 'seismic'
% data_type: 'radar'
% data_type: 'cluster'

eval(pfile);       % run script with the input parameters
addpath('seismic/');
if(numel(it)==1)
    modelvec = load(sprintf('%s%i.txt',d.ds(1).name,it));
else
    modelvec = it;
end
mod = reshape(modelvec,d.grid.my,d.grid.mx,d.grid.mz);
mod(:,:,end+1) = 0;%mod(:,:,end);
mod(:,end+1,:) = 0;%mod(:,end,:);
mod(end+1,:,:) = 0;%mod(end,:,:);


% xshift = 6.65;
% yshift = 75.75;

xshift = 6.65;
yshift = 75.75;

x = (0:d.grid.dx:d.grid.dx*d.grid.mx) + xshift;
y = (0:d.grid.dy:d.grid.dy*d.grid.my) + yshift;
z = 0:d.grid.dz:d.grid.dz*d.grid.mz;

%figure; hold on
V = slice(x,y,z,mod(:,:,:),xp,yp,zp);
% mod(x<45,:,:) = NaN;
% mod(:,y>120,:) = NaN;
% mod(:,y<100,:) = NaN;
if(nargin>5 && iso)
    p = patch(isosurface(x,y,z,mod,-30),'facecolor','y');
end
set(V,'edgecolor','none')

% set(gca,'fontsize',20,'linewidth',2)
set(gca,'dataaspectratio',[1 1 1],'layer','top')
%    set(gca,'xtick',[0.92 2.92 4.92],'xticklabel',[0 2 4])
%    set(gca,'ytick',[1.7 3.7 5.7],'yticklabel',[4 2 0])
%    set(gca,'ztick',0:2:10,'zticklabel',0:2:10,'zdir','reverse')
% set(gca,'xdir','reverse');
xlabel('Northing [m]');
ylabel('Easting [m]')
zlabel('Elevation [m]');
axis equal tight
    
Vp = mod;
x_seis = x;
y_seis = y;
z_seis = z;
%caxis([4.7 5.5]);

% 
% scoords = load(sprintf('%s/stations.dat',d.ds(1).folder));
% rcoords = load(sprintf('%s/shots.dat',d.ds(1).folder));
% 
% plot3(scoords(:,3)+yshift,scoords(:,2)+xshift,scoords(:,4),'r+');
% plot3(rcoords(:,3)+yshift,rcoords(:,2)+xshift,rcoords(:,4),'ko');

grid on
box on
axis equal tight
view(90-25,30)
%view(209,46)
xlim;
ylim;
zlim;

