%% C63 self-contained two-body demonstration
% Reproduces the spacecraft maneuver geometry, not debris avoidance results.
% No TLE file needed. MATLAB R2026b; base MATLAB ODE solver.
clear; clc; close all;
mu = 3.986004418e14; Re = 6378137;
altitude = 500e3; inclination = deg2rad(51.6);
r0 = [Re+altitude;0;0];
vCircular = sqrt(mu/norm(r0));
v0 = [0;vCircular*cos(inclination);vCircular*sin(inclination)];
y0 = [r0;v0];
t1 = 2600; t2 = 8000; tEnd = 10800;
dv1 = -0.40; dv2 = 0.20;
opts = odeset('RelTol',1e-11,'AbsTol',1e-7);
s0 = ode113(@(t,y)rhs(t,y,mu),[0 tEnd],y0,opts);
s1 = ode113(@(t,y)rhs(t,y,mu),[0 t1],y0,opts);
y1 = deval(s1,t1); y1(4:6) = y1(4:6)+dv1*tangent(y1);
s2 = ode113(@(t,y)rhs(t,y,mu),[t1 t2],y1,opts);
y2 = deval(s2,t2); y2(4:6) = y2(4:6)+dv2*tangent(y2);
s3 = ode113(@(t,y)rhs(t,y,mu),[t2 tEnd],y2,opts);
yFinal = deval(s3,tEnd);
a0 = semiMajorAxis(y0,mu);
aFinal = semiMajorAxis(yFinal,mu);
eFinal = eccentricity(yFinal,mu);
fprintf('C63 total DV: %.3f m/s\n',abs(dv1)+abs(dv2));
fprintf('Final |Delta a|: %.6f km\n',abs(aFinal-a0)/1e3);
fprintf('Final eccentricity: %.9f\n',eFinal);
fprintf('Orbit tolerance pass: %d\n',abs(aFinal-a0)<=5e3 && eFinal<=1e-3);
tt = linspace(0,tEnd,2000);
p0 = deval(s0,tt,1:3)/1e3;
p = zeros(3,numel(tt));
i1=tt<t1; i2=tt>=t1 & tt<t2; i3=tt>=t2;
p(:,i1)=deval(s1,tt(i1),1:3)/1e3;
p(:,i2)=deval(s2,tt(i2),1:3)/1e3;
p(:,i3)=deval(s3,tt(i3),1:3)/1e3;
figure('Color','w'); hold on; grid on; axis equal;
plot3(p0(1,:),p0(2,:),p0(3,:),'--','LineWidth',1);
plot3(p(1,:),p(2,:),p(3,:),'LineWidth',1.5);
b1=deval(s1,t1,1:3)/1e3; b2=deval(s2,t2,1:3)/1e3;
scatter3(b1(1),b1(2),b1(3),70,'filled');
scatter3(b2(1),b2(2),b2(3),70,'filled');
xlabel('Inertial x (km)'); ylabel('Inertial y (km)'); zlabel('Inertial z (km)');
title('C63: two tangential burns (spacecraft only)');
legend('No maneuver','C63','Burn 1: -0.40 m/s','Burn 2: +0.20 m/s','Location','best');
view(35,22);
% This demo cannot verify closest approach to any debris object.
function dy=rhs(~,y,mu)
    r=y(1:3); dy=[y(4:6);-mu*r/norm(r)^3];
end
function T=tangent(y)
    r=y(1:3); v=y(4:6);
    R=r/norm(r); N=cross(r,v); N=N/norm(N);
    T=cross(N,R); T=T/norm(T);
end
function a=semiMajorAxis(y,mu)
    r=y(1:3);v=y(4:6);
    a=1/(2/norm(r)-dot(v,v)/mu);
end
function e=eccentricity(y,mu)
    r=y(1:3);v=y(4:6);
    e=norm(((dot(v,v)-mu/norm(r))*r-dot(r,v)*v)/mu);
end
