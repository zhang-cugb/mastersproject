function [SZ_tunnel] = ThreeD_Shear_zones_tunnels(SZ_tunnel,iv)
% function [ output_args ] = ThreeD_Fracture_discs( OPTVLogs )
figure(iv)
Tunnel_disc(:,1:3) = [1.*sind(SZ_tunnel.Azimuth+180).*sind(SZ_tunnel.Dip), 1.*cosd(SZ_tunnel.Azimuth+180).*sind(SZ_tunnel.Dip),    1.*cosd(SZ_tunnel.Azimuth)];
for ii = 1:length(SZ_tunnel.Tunnel)
    if SZ_tunnel.SZ(ii) == 11 || SZ_tunnel.SZ(ii) == 12 || SZ_tunnel.SZ(ii) == 13
        c = 'r-';
    elseif SZ_tunnel.SZ(ii) == 31 || SZ_tunnel.SZ(ii) == 32
        c = 'g-';
    end
        radius= 2;
        center= [SZ_tunnel.x(ii) SZ_tunnel.y(ii), SZ_tunnel.z(ii)];
        normal = [Tunnel_disc(ii,1), Tunnel_disc(ii,2), Tunnel_disc(ii,3)];
        theta=0:0.01:2*pi;
        v=null(normal);
        points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
        patch(points(1,:),points(2,:),points(3,:),c);
        
end

end

