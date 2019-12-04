clear
close all

% Import Image
Image = importdata('AU_Mala41.png');

% Create new variables in the base workspace from those fields.
vars = fieldnames(Image);
for i = 1:length(vars)
    assignin('base', vars{i}, Image.(vars{i}));
end

% extract colorinformation
c1 = cdata(:,:,1);
c2 = cdata(:,:,2);
c3 = cdata(:,:,3);

% Estimate length of picture
sx = 1:length(cdata(1,:,1));
sy = 1:length(cdata(:,1,1));
lx = max(sx);
ly = max(sy);

% Length of picture in [m]
x = 40;
y = 20;

% Dip angle of image
Dip = 45;

% Pixel size
dx = x/lx;
dy = y/ly;

rx = dx:dx:x;
ry = dy:dy:y;
z = ones(length(ry),length(rx));

for ii = 1: length(z(1,:))
    z(:,ii) = z(:,ii)*rx(ii)*sind(Dip);
    x_final(ii) = rx(ii)*cosd(Dip);
end
y_final = ry;

figure
imagesc(rx,ry,cdata)
figure
s=surf(rx,ry,z)
s.CData = cdata(:,:,[1 2 3])
set(s,'edgecolor','none')

figure
s=surf(x_final,y_final,z)
s.CData = cdata(:,:,[1 2 3])
set(s,'edgecolor','none')