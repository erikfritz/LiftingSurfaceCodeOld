function [R] = plotGeometry(d1xLE,d1yLE,d1zLE,d1xTE,d1yTE,d1zTE,d1CPgridX,d1CPgridY,d1CPgridZ,...
    d2gridX,d2gridY,d2gridZ,d2VRgridX,d2VRgridY,d2VRgridZ,d2wakeX,d2wakeY,d2wakeZ,d1Gamma,velX,dt,Nt)

d2Gamma = reshape(d1Gamma,size(d2wakeX,1)-1,size(d2wakeX,2)-1);
% figure('Name','Geometry','NumberTitle','off')
figure(1)
clf
set(gcf, 'Position', get(0, 'Screensize'));
hold on; grid on;
plot3(d1xLE,d1yLE,d1zLE,'k')
plot3(d1xTE,d1yTE,d1zTE,'k')
scatter3(d1CPgridX,d1CPgridY,d1CPgridZ,'xk')
surf(d2gridX,d2gridY,d2gridZ,'FaceColor','w')
surf(d2VRgridX,d2VRgridY,d2VRgridZ,'FaceColor','none','EdgeColor','b')
surf(d2wakeX,d2wakeY,d2wakeZ,d2Gamma,'FaceAlpha',0.5)

xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
xlim([-1 velX*dt*Nt+2])
c = colorbar('southoutside');
colormap jet
c.Label.String = 'Circulation \Gamma [m^2/s]';
view([5760 80])
daspect([1 1 1])
drawnow

pause(0.05)
R = [];
% print(R)
end