function plot_GPR(ivGPR)

Image = importdata('GPR_AU.png');

% Create new variables in the base workspace from those fields.
% vars = fieldnames(Image);
% for i = 1:length(vars)
%     assignin('base', vars{i}, Image.(vars{i}));
% end
cdata = Image;
% extract colorinformation
c1 = cdata(:,:,1);
c2 = cdata(:,:,2);
c3 = cdata(:,:,3);

% Estimate length of picture
sx = 1:length(cdata(:,1,1));
sy = 1:length(cdata(1,:,1));
lx = max(sx);
ly = max(sy);

% Length of picture in [m]
x = 40;
y = -52;

% Dip angle of image
Dip = -60;

% Pixel size
dx = x/lx;
dy = y/ly;

rx = dx:dx:x;
ry = dy:dy:y;
z = ones(length(ry),length(rx));

for ii = 1: length(z(1,:))
    z(:,ii) = z(:,ii)*ry(ii)*sind(Dip);
    x_final(ii) = rx(ii)*cosd(Dip);
end
y_final = ry;

figure(ivGPR)
s=surf(x_final+55,y_final+145,z+34-24.5);
s.CData = rot90(cdata,3);
set(s,'edgecolor','none');
material dull
axis equal
xlabel('east')
end

