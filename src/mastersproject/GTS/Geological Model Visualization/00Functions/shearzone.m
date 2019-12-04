function S12_BH = shearzone(file, pathSH, no, BH, GTS_coordinates)
% CALCULATE DIP AND DIP DIRCTOIN FROM 3 GIVEN POINT AT AU & VE TUNNEL
shearcircle     = load_S1S3(file, pathSH);
intersection    = loadtunnelintersection('Tunnel_intersections.txt', pathSH);

for i=1:size(shearcircle.BH,1)
    S12_BH.depth(i) = shearcircle.depth(i);
    S12_BH.az(i)    = shearcircle.az(i);
    S12_BH.dip(i)   = shearcircle.dip(i);
    S12_BH.ap(i)    = shearcircle.ap(i);
    BHmouth = [BH.coordinates(i).x_m; BH.coordinates(i).y_m; BH.coordinates(i).z_m];
    BHend   = [BH.coordinates(i).x_e; BH.coordinates(i).y_e; BH.coordinates(i).z_e];
    n1 = (BHend-BHmouth); n1=n1/norm(n1);
    S12_BH.x_m(i) = BHmouth(1)+n1(1)*shearcircle.depth(i);
    S12_BH.y_m(i) = BHmouth(2)+n1(2)*shearcircle.depth(i);
    S12_BH.z_m(i) = BHmouth(3)+n1(3)*shearcircle.depth(i);      
end
for ii = 1:size(intersection.she,1)
    if intersection.she(ii) == no
        i = i+1;
        S12_BH.depth(i) = NaN;
        S12_BH.az(i)    = intersection.az(ii);
        S12_BH.dip(i)   = intersection.dip(ii);
        S12_BH.ap(i)    = NaN;
        S12_BH.x_m(i) = intersection.x(ii)-GTS_coordinates.x;
        S12_BH.y_m(i) = intersection.y(ii)-GTS_coordinates.y;
        S12_BH.z_m(i) = intersection.z(ii)-GTS_coordinates.z;
    end
end
S12_BH.depth= S12_BH.depth';
S12_BH.az   = S12_BH.az';
S12_BH.dip  = S12_BH.dip';
S12_BH.ap   = S12_BH.ap';
S12_BH.x_m  = S12_BH.x_m';
S12_BH.y_m  = S12_BH.y_m';
S12_BH.z_m  = S12_BH.z_m';

end