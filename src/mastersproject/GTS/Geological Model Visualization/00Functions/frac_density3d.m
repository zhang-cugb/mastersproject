function [ FracturePerMeter ] = frac_density3d( FBS, SBH, INJ, GEO, PRP ,iv )

% Import Fracture per Meter log
[ FracturePerMeter ] = importFracDensity;

count = 1;
for ii = 1 : length(FracturePerMeter.Depthm(:,1))
    if FracturePerMeter.Borehole(ii,1) == 'INJ1'
        FracturePerMeter.x(count,1) = INJ(1,1)+FracturePerMeter.Depthm(ii).*sind(INJ(1,6)).*cosd(INJ(1,7));
        FracturePerMeter.y(count,1) = INJ(1,2)+FracturePerMeter.Depthm(ii).*cosd(INJ(1,6))*cosd(INJ(1,7));
        FracturePerMeter.z(count,1) = INJ(1,3)+FracturePerMeter.Depthm(ii).*sind(INJ(1,7));
    elseif FracturePerMeter.Borehole(ii,1) == 'INJ2'
        FracturePerMeter.x(count,1) = [INJ(2,1)+FracturePerMeter.Depthm(ii).*sind(INJ(2,6)).*cosd(INJ(2,7))];
        FracturePerMeter.y(count,1) = [INJ(2,2)+FracturePerMeter.Depthm(ii).*cosd(INJ(2,6))*cosd(INJ(2,7))];
        FracturePerMeter.z(count,1) = [INJ(2,3)+FracturePerMeter.Depthm(ii).*sind(INJ(2,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'FBS1'
        FracturePerMeter.x(count,1) = [FBS(1,1)+FracturePerMeter.Depthm(ii).*sind(FBS(1,6)).*cosd(FBS(1,7))];
        FracturePerMeter.y(count,1) = [FBS(1,2)+FracturePerMeter.Depthm(ii).*cosd(FBS(1,6))*cosd(FBS(1,7))];
        FracturePerMeter.z(count,1) = [FBS(1,3)+FracturePerMeter.Depthm(ii).*sind(FBS(1,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'FBS2'
        FracturePerMeter.x(count,1) = [FBS(2,1)+FracturePerMeter.Depthm(ii).*sind(FBS(2,6)).*cosd(FBS(2,7))];
        FracturePerMeter.y(count,1) = [FBS(2,2)+FracturePerMeter.Depthm(ii).*cosd(FBS(2,6))*cosd(FBS(2,7))];
        FracturePerMeter.z(count,1) = [FBS(2,3)+FracturePerMeter.Depthm(ii).*sind(FBS(2,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'FBS3'
        FracturePerMeter.x(count,1) = [FBS(3,1)+FracturePerMeter.Depthm(ii).*sind(FBS(3,6)).*cosd(FBS(3,7))];
        FracturePerMeter.y(count,1) = [FBS(3,2)+FracturePerMeter.Depthm(ii).*cosd(FBS(3,6))*cosd(FBS(3,7))];
        FracturePerMeter.z(count,1) = [FBS(3,3)+FracturePerMeter.Depthm(ii).*sind(FBS(3,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'PRP1'
        FracturePerMeter.x(count,1) = [PRP(1,1)+FracturePerMeter.Depthm(ii).*sind(PRP(1,6)).*cosd(PRP(1,7))];
        FracturePerMeter.y(count,1) = [PRP(1,2)+FracturePerMeter.Depthm(ii).*cosd(PRP(1,6))*cosd(PRP(1,7))];
        FracturePerMeter.z(count,1) = [PRP(1,3)+FracturePerMeter.Depthm(ii).*sind(PRP(1,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'PRP2'
        FracturePerMeter.x(count,1) = [PRP(2,1)+FracturePerMeter.Depthm(ii).*sind(PRP(2,6)).*cosd(PRP(2,7))];
        FracturePerMeter.y(count,1) = [PRP(2,2)+FracturePerMeter.Depthm(ii).*cosd(PRP(2,6))*cosd(PRP(2,7))];
        FracturePerMeter.z(count,1) = [PRP(2,3)+FracturePerMeter.Depthm(ii).*sind(PRP(2,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'PRP3'
        FracturePerMeter.x(count,1) = [PRP(3,1)+FracturePerMeter.Depthm(ii).*sind(PRP(3,6)).*cosd(PRP(3,7))];
        FracturePerMeter.y(count,1) = [PRP(3,2)+FracturePerMeter.Depthm(ii).*cosd(PRP(3,6))*cosd(PRP(3,7))];
        FracturePerMeter.z(count,1) = [PRP(3,3)+FracturePerMeter.Depthm(ii).*sind(PRP(3,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'SBH1'
        FracturePerMeter.x(count,1) = [SBH(1,1)+FracturePerMeter.Depthm(ii).*sind(SBH(1,6)).*cosd(SBH(1,7))];
        FracturePerMeter.y(count,1) = [SBH(1,2)+FracturePerMeter.Depthm(ii).*cosd(SBH(1,6))*cosd(SBH(1,7))];
        FracturePerMeter.z(count,1) = [SBH(1,3)+FracturePerMeter.Depthm(ii).*sind(SBH(1,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'SBH3'
        FracturePerMeter.x(count,1) = [SBH(2,1)+FracturePerMeter.Depthm(ii).*sind(SBH(2,6)).*cosd(SBH(2,7))];
        FracturePerMeter.y(count,1) = [SBH(2,2)+FracturePerMeter.Depthm(ii).*cosd(SBH(2,6))*cosd(SBH(2,7))];
        FracturePerMeter.z(count,1) = [SBH(2,3)+FracturePerMeter.Depthm(ii).*sind(SBH(2,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'SBH4'
        FracturePerMeter.x(count,1) = [SBH(3,1)+FracturePerMeter.Depthm(ii).*sind(SBH(3,6)).*cosd(SBH(3,7))];
        FracturePerMeter.y(count,1) = [SBH(3,2)+FracturePerMeter.Depthm(ii).*cosd(SBH(3,6))*cosd(SBH(3,7))];
        FracturePerMeter.z(count,1) = [SBH(3,3)+FracturePerMeter.Depthm(ii).*sind(SBH(3,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'GEO1'
        FracturePerMeter.x(count,1) = [GEO(1,1)+FracturePerMeter.Depthm(ii).*sind(GEO(1,6)).*cosd(GEO(1,7))];
        FracturePerMeter.y(count,1) = [GEO(1,2)+FracturePerMeter.Depthm(ii).*cosd(GEO(1,6))*cosd(GEO(1,7))];
        FracturePerMeter.z(count,1) = [GEO(1,3)+FracturePerMeter.Depthm(ii).*sind(GEO(1,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'GEO2'
        FracturePerMeter.x(count,1) = [GEO(2,1)+FracturePerMeter.Depthm(ii).*sind(GEO(2,6)).*cosd(GEO(2,7))];
        FracturePerMeter.y(count,1) = [GEO(2,2)+FracturePerMeter.Depthm(ii).*cosd(GEO(2,6))*cosd(GEO(2,7))];
        FracturePerMeter.z(count,1) = [GEO(2,3)+FracturePerMeter.Depthm(ii).*sind(GEO(2,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'GEO3'
        FracturePerMeter.x(count,1) = [GEO(3,1)+FracturePerMeter.Depthm(ii).*sind(GEO(3,6)).*cosd(GEO(3,7))];
        FracturePerMeter.y(count,1) = [GEO(3,2)+FracturePerMeter.Depthm(ii).*cosd(GEO(3,6))*cosd(GEO(3,7))];
        FracturePerMeter.z(count,1) = [GEO(3,3)+FracturePerMeter.Depthm(ii).*sind(GEO(3,7))];
    elseif FracturePerMeter.Borehole(ii,1) == 'GEO4'
        FracturePerMeter.x(count,1) = [GEO(4,1)+FracturePerMeter.Depthm(ii).*sind(GEO(4,6)).*cosd(GEO(4,7))];
        FracturePerMeter.y(count,1) = [GEO(4,2)+FracturePerMeter.Depthm(ii).*cosd(GEO(4,6))*cosd(GEO(4,7))];
        FracturePerMeter.z(count,1) = [GEO(4,3)+FracturePerMeter.Depthm(ii).*sind(GEO(4,7))];
    end
    count = count+1;
end


figure(iv)
scatter3(FracturePerMeter.x(:,1),FracturePerMeter.y(:,1),FracturePerMeter.z(:,1),40,FracturePerMeter.FracCount(:,1),'filled');
caxis([0 15]);
colormap('hot')
% 
% figure(100)
% plot(FracturePerMeter.Vp,FracturePerMeter.FracCount ,'*')
% xlabel('Vp [km/s]')
% ylabel('Frac/meter [1/m]')
% % title('Fracture Density')
% 
% figure(101)
% plot(FracturePerMeter.Vp_averaged(1:end-2),values(:,4) ,'*')
% xlabel('Vp [km/s]')
% ylabel('averaged Frac/meter [1/m]')
% % title('Averaged Fracture Density')
end

