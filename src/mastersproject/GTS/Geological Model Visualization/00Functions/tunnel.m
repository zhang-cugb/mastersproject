function tunnel( GTS_coordinates, num_fig)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%Transparency of tunnels
Alpha = 1;
%% PATH ABSOLUTE TO LOADED DATA
% pathT = ('..\01 GTS\02_TunnelCavernCoordinates\');
% LOAD Tunnel Points  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tunnel.VE = load_tunnel('VETunnel.txt',[]);
Tunnel.AU = load_tunnel('AUTunnel.txt',[]);
Tunnel.AUcavern = excav_points('AUcavern.txt',[], 2);
Tunnel.AUexcav = excav_points('AUexcav.txt', [],2);
Tunnel.VEcavern = excav_points('VEcavern.txt',[], 2);
Tunnel.AUgallery = excav_points('AUgallery.txt',[], 6);

Tunnel.VE.diameter = 3.495; %m
Tunnel.AU.diameter = 3.475; %m
Tunnel.AUgallery.diameter = 2.4; %m
Tunnel.AUgallery.lengthNS = 67.0; %m
Tunnel.AUgallery.lengthEW= 6.7; %m
Tunnel.AUgallery.distfromStoCross = 32.33; %m

Tunnel.AU.x_m   = Tunnel.AU.x_m - GTS_coordinates.x;
Tunnel.AU.y_m   = Tunnel.AU.y_m - GTS_coordinates.y;
Tunnel.AU.z_m   = Tunnel.AU.z_m - GTS_coordinates.z;
Tunnel.VE.x_m   = Tunnel.VE.x_m - GTS_coordinates.x;
Tunnel.VE.y_m   = Tunnel.VE.y_m - GTS_coordinates.y;
Tunnel.VE.z_m   = Tunnel.VE.z_m - GTS_coordinates.z;
Tunnel.AUgallery.x_m   = Tunnel.AUgallery.x_m - GTS_coordinates.x;
Tunnel.AUgallery.y_m   = Tunnel.AUgallery.y_m - GTS_coordinates.y;
Tunnel.AUgallery.z_m   = Tunnel.AUgallery.z_m - GTS_coordinates.z;
Tunnel.AUcavern.x_m   = Tunnel.AUcavern.x_m - GTS_coordinates.x;
Tunnel.AUcavern.y_m   = Tunnel.AUcavern.y_m - GTS_coordinates.y;
Tunnel.AUcavern.z_m   = Tunnel.AUcavern.z_m - GTS_coordinates.z;
Tunnel.VEcavern.x_m   = Tunnel.VEcavern.x_m - GTS_coordinates.x;
Tunnel.VEcavern.y_m   = Tunnel.VEcavern.y_m - GTS_coordinates.y;
Tunnel.VEcavern.z_m   = Tunnel.VEcavern.z_m - GTS_coordinates.z;
Tunnel.AUexcav.x_m   = Tunnel.AUexcav.x_m - GTS_coordinates.x;
Tunnel.AUexcav.y_m   = Tunnel.AUexcav.y_m - GTS_coordinates.y;
Tunnel.AUexcav.z_m   = Tunnel.AUexcav.z_m - GTS_coordinates.z;


for ii = 1:num_fig
    
    figure(ii)
% PLOT TUNNEL
% construct first point for VE tunnel 
    i=1;
    P1 = [Tunnel.VE.x_m(i); Tunnel.VE.y_m(i); Tunnel.VE.z_m(i)];
    i=2;
    P2 = [Tunnel.VE.x_m(i); Tunnel.VE.y_m(i); Tunnel.VE.z_m(i)];
    n = (P2-P1); n=n/norm(n);
    %norm([Tunnel.VE.x_m(3)-Tunnel.VE.x_m(2);Tunnel.VE.y_m(3)-Tunnel.VE.y_m(2);Tunnel.VE.z_m(3)-Tunnel.VE.z_m(2)])
    Tunnel.VE.x_m = [Tunnel.VE.x_m; Tunnel.VE.x_m(2)+n(1)*11/2];
    Tunnel.VE.y_m = [Tunnel.VE.y_m; Tunnel.VE.y_m(2)+n(2)*11/2];
    Tunnel.VE.z_m = [Tunnel.VE.z_m; Tunnel.VE.z_m(2)+n(3)*11/2];

    % constructs VE tunnel
    i=1;
    Tunnel.VE.Cylinder(i) = Cylinder([Tunnel.VE.x_m(i); Tunnel.VE.y_m(i);Tunnel.VE.z_m(i)], ...
        [Tunnel.VE.x_m(end); Tunnel.VE.y_m(end);Tunnel.VE.z_m(end)], Tunnel.VE.diameter/2, ...
        20, [.3 1 1], 1,0);

    % constructs VE cavern
    %scatter3(Tunnel.VEcavern.x_m, Tunnel.VEcavern.y_m, Tunnel.VEcavern.z_m);
    hold on
    v = [Tunnel.VEcavern.x_m, Tunnel.VEcavern.y_m, Tunnel.VEcavern.z_m];
    Tunnel.VE.face = [9 8 11; 11 10 8;13 15 9; 9 11 13; 15 13 12; 12 15 14;...
        5 8 9; 5 8 10; 4 10 5; 4 10 11; 11 13 3; 4 3 11;...
        7 12 13; 13 7 3; 1 4 3; 4 1 5; 2 9 15; 9 5 2; 1 2 5;...
        6 14 15; 6 2 15; 1 3 2; 3 7 6; 6 2 3; 14 12 6; 7 12 6];
    patch('Faces',Tunnel.VE.face,'Vertices',v,'FaceColor',[.3 1 1],'FaceAlpha', 0.5 , 'EdgeColor', 'none')
    %Tunnel.VEcavern.TRI = delaunay(Tunnel.VEcavern.x_m, Tunnel.VEcavern.y_m, Tunnel.VEcavern.z_m);
    %Tunnel.VEcavern.Surf= trisurf(Tunnel.VEcavern.TRI,Tunnel.VEcavern.x_m, Tunnel.VEcavern.y_m, Tunnel.VEcavern.z_m, 'FaceColor', [.3 1 1],'FaceAlpha', 0.5);
    
    % constructs last point for AU tunnel
    i=1;
    P1 = [Tunnel.AU.x_m(i); Tunnel.AU.y_m(i); Tunnel.AU.z_m(i)];
    i=2;
    P2 = [Tunnel.AU.x_m(i); Tunnel.AU.y_m(i); Tunnel.AU.z_m(i)];
    n = (P2-P1); n=n/norm(n);
    Tunnel.AU.x_m = [Tunnel.AU.x_m; Tunnel.AU.x_m(i)+n(1)*(Tunnel.AU.length(i)/2+Tunnel.AU.extra_l(i))-.5];
    Tunnel.AU.y_m = [Tunnel.AU.y_m; Tunnel.AU.y_m(i)+n(2)*(Tunnel.AU.length(i)/2+Tunnel.AU.extra_l(i))];
    Tunnel.AU.z_m = [Tunnel.AU.z_m; Tunnel.AU.z_m(i)+n(3)*(Tunnel.AU.length(i)/2+Tunnel.AU.extra_l(i))];
    % constructs AU tunnel
    i=1;
    Tunnel.AU.Cylinder(i) = Cylinder([Tunnel.AU.x_m(i); Tunnel.AU.y_m(i);Tunnel.AU.z_m(i)], ...
        [Tunnel.AU.x_m(end); Tunnel.AU.y_m(end);Tunnel.AU.z_m(end)], Tunnel.AU.diameter/2, ...
        20, [.3 1 1], 1,0);
    % constructs AU cavern
    %scatter3(Tunnel.AUcavern.x_m, Tunnel.AUcavern.y_m, Tunnel.AUcavern.z_m)
    hold on
    v = [Tunnel.AUcavern.x_m, Tunnel.AUcavern.y_m, Tunnel.AUcavern.z_m];
    Tunnel.AUcavern.face = [22 23 21; 21 20 24; 23 24 21;...
        29 11 20; 11 10 20; 20 10 24; 18 10 24; 26 11 10; 10 18 25; 26 25 10; 26 29 11;...
        20 21 29; 21 29 30; 21 22 30; 28 22 30; 28 27 30; 2 30 29; 29 26 2; 26 2 25; 2 27 30;...
        16 12 23; 23 12 22; 12 28 22; 8 12 16; 12 27 28; 12 6 27; 2 27 6; 2 6 25;...
        18 9 24; 9 23 24; 9 16 23; 9 8 16; 25 18 9; 25 9 6; 6 9 8; 6 8 12];
    patch('Faces',Tunnel.AUcavern.face,'Vertices',v,'FaceColor',[.3 1 1],'FaceAlpha', 0.5 , 'EdgeColor', 'none')
       
    %Tunnel.AUcavern.TRI = delaunay(Tunnel.AUcavern.x_m, Tunnel.AUcavern.y_m, Tunnel.AUcavern.z_m);
    %Tunnel.AUcavern.Surf= trisurf(Tunnel.AUcavern.TRI,Tunnel.AUcavern.x_m, Tunnel.AUcavern.y_m, Tunnel.AUcavern.z_m, 'FaceColor', [.3 1 1],'FaceAlpha', 0.5);
    
    
    % constructs AU excavation
    %scatter3(Tunnel.AUexcav.x_m, Tunnel.AUexcav.y_m, Tunnel.AUexcav.z_m)
    hold on
    v = [Tunnel.AUexcav.x_m, Tunnel.AUexcav.y_m, Tunnel.AUexcav.z_m];
    Tunnel.AUexcav.face = [7 9 10; 6 7 10; 16 6 10; 10 11 16; 15 11 16; 6 16 15;...
        8 7 9; 8 12 9; 4 8 7; 14 6 15; 3 12 8; 3 8 4; 3 2 4; 13 4 2;...
        5 2 1; 13 5 2; 13 14 4; 5 1 11; 5 13 11; 13 11 15; 15 14 13;...
        4 7 6; 4 6 14; 12 9 10; 11 12 10; 11 8 12; 1 11 2; 11 2 8; 2 4 8];
    patch('Faces',Tunnel.AUexcav.face,'Vertices',v,'FaceColor',[.3 1 1],'FaceAlpha', 0.5 , 'EdgeColor', 'none')
    %Tunnel.AUexcav.TRI = delaunay(Tunnel.AUexcav.x_m, Tunnel.AUexcav.y_m, Tunnel.AUexcav.z_m);
    %Tunnel.AUexcav.Surf= trisurf(Tunnel.AUexcav.TRI,Tunnel.AUexcav.x_m, Tunnel.AUexcav.y_m, Tunnel.AUexcav.z_m, 'FaceColor', [.3 1 1],'FaceAlpha', 0.5);


    % constructs AU gallery
    Tunnel.AUgallery.Cylinder(1) = Cylinder([Tunnel.AUgallery.x_m(1); Tunnel.AUgallery.y_m(1);Tunnel.AUgallery.z_m(1)], ...
        [Tunnel.AUgallery.x_m(2); Tunnel.AUgallery.y_m(2);Tunnel.AUgallery.z_m(2)], Tunnel.AUgallery.diameter/2, ...
        20, [.3 1 1], 1,0);
    i=1;
    P1 = [Tunnel.AUgallery.x_m(i); Tunnel.AUgallery.y_m(i); Tunnel.AUgallery.z_m(i)];
    i=2;
    P2 = [Tunnel.AUgallery.x_m(i); Tunnel.AUgallery.y_m(i); Tunnel.AUgallery.z_m(i)];
    n = (P2-P1); n=n/norm(n);
    i=3;
    P3 = [Tunnel.AUgallery.x_m(i); Tunnel.AUgallery.y_m(i); Tunnel.AUgallery.z_m(i)];
    i=4;
    P4 = [Tunnel.AUgallery.x_m(1)+n(1)*Tunnel.AUgallery.distfromStoCross; Tunnel.AUgallery.y_m(1)+n(2)*Tunnel.AUgallery.distfromStoCross; Tunnel.AUgallery.z_m(1)+n(3)*Tunnel.AUgallery.distfromStoCross];
    n = (P4-P3); n=n/norm(n);

    Tunnel.AUgallery.x_m = [Tunnel.AUgallery.x_m; Tunnel.AUgallery.x_m(3)+n(1)*Tunnel.AUgallery.lengthEW(1)/2];
    Tunnel.AUgallery.y_m = [Tunnel.AUgallery.y_m; Tunnel.AUgallery.y_m(3)+n(2)*Tunnel.AUgallery.lengthEW(1)/2];
    Tunnel.AUgallery.z_m = [Tunnel.AUgallery.z_m; Tunnel.AUgallery.z_m(3)+n(3)*Tunnel.AUgallery.lengthEW(1)/2];
    Tunnel.AUgallery.Cylinder(2) = Cylinder([Tunnel.AUgallery.x_m(3); Tunnel.AUgallery.y_m(3);Tunnel.AUgallery.z_m(3)], ...
        [Tunnel.AUgallery.x_m(4); Tunnel.AUgallery.y_m(4);Tunnel.AUgallery.z_m(4)], Tunnel.AUgallery.diameter/2, ...
        20, [.3 1 1], 1,0);
    % plot3(Tunnel.VE.x_m, Tunnel.VE.y_m, Tunnel.VE.z_m, '-k')
    % plot3(Tunnel.AU.x_m, Tunnel.AU.y_m, Tunnel.AU.z_m, '-k')
    % scatter3(Tunnel.VE.x_m, Tunnel.VE.y_m, Tunnel.VE.z_m, 'xk')
    % scatter3(Tunnel.AU.x_m, Tunnel.AU.y_m, Tunnel.AU.z_m, 'xk')

material shiny
end
end

