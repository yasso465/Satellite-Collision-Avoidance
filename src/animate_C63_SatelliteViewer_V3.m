%% C63 V3 — native MATLAB satelliteScenarioViewer
% MATLAB R2026b + Aerospace Toolbox. Run from repository ROOT.
% Reads actual TLEs; no synthetic debris trajectories.
% Quick preview: set exportMP4=false. Rendering may take a few minutes.
clear; clc; close all;
frameStep = 10; % seconds per scenario step
noradIDs = ["44804","54771"]; % critical retained-catalog encounters
startUTC = datetime(2026,10,12,0,0,0,'TimeZone','UTC');
durationS = 10800; mu = 3.986004418e14; Re = 6378137;
% Resolve paths relative to this script, not MATLAB's current directory.
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
tlePath = fullfile(projectRoot,'data','debris_matlab.tle');
assert(isfile(tlePath), 'TLE missing at: %s', tlePath);
% Extract only the requested actual TLE records, preserving 3-line format.
lines = readlines(tlePath); lines = lines(strlength(strtrim(lines))>0);
assert(mod(numel(lines),3)==0,'Expected 3-line TLE catalog.');
tmpFiles = strings(1,numel(noradIDs));
for k=1:numel(noradIDs)
    found=false;
    for j=1:3:numel(lines)
        l1=char(lines(j+1));
        if numel(l1)>=7 && string(strtrim(l1(3:7)))==noradIDs(k)
            tmpFiles(k)=string([tempname '.tle']);
            fid=fopen(tmpFiles(k),'w');
            assert(fid>0,'Unable to create temporary TLE');
            fprintf(fid,'%s\n%s\n%s\n',char(lines(j)),char(lines(j+1)),char(lines(j+2)));
            fclose(fid); found=true; break;
        end
    end
    assert(found,'NORAD %s not found in %s',noradIDs(k),tlePath);
end
cleanupObj=onCleanup(@()deleteTemps(tmpFiles));
% Exact spacecraft initial state follows the validation workflow.
sc0=satelliteScenario(startUTC,startUTC+seconds(durationS),frameStep);
mission=satellite(sc0,Re+500e3,0,51.6,0,0,0,'OrbitPropagator','two-body-keplerian');
[p0,v0]=states(mission,'CoordinateFrame','inertial');
y0=[p0(:,1);v0(:,1)];
t1=2600; t2=8000; dv1=-0.40; dv2=+0.20;
opts=odeset('RelTol',1e-10,'AbsTol',1e-5);
s1=ode113(@(~,y)rhs(y,mu),[0 t1],y0,opts);
y1=deval(s1,t1);y1(4:6)=y1(4:6)+dv1*tangent(y1);
s2=ode113(@(~,y)rhs(y,mu),[t1 t2],y1,opts);
y2=deval(s2,t2);y2(4:6)=y2(4:6)+dv2*tangent(y2);
s3=ode113(@(~,y)rhs(y,mu),[t2 durationS],y2,opts);
sols={s1,s2,s3};
% One actual SGP4 propagation for each selected object.
t=0:frameStep:durationS; N=numel(t);
satP=positionsAt(sols,t,t1,t2);
debP=nan(3,N,numel(noradIDs));
for k=1:numel(noradIDs)
    sc=satelliteScenario(startUTC,startUTC+seconds(durationS),frameStep);
    d=satellite(sc,char(tmpFiles(k)),'OrbitPropagator','sgp4');
    [p,~,ts]=states(d,'CoordinateFrame','inertial');
    dt=seconds(ts(:)-startUTC)';
    assert(size(p,2)==N && max(abs(dt-t))<1e-5, ...
        'SGP4 sample grid mismatch for NORAD %s',noradIDs(k));
    debP(:,:,k)=p;
end
distKm=zeros(numel(noradIDs),N);
for k=1:numel(noradIDs)
    distKm(k,:)=vecnorm(satP-debP(:,:,k),2,1)/1000;
    fprintf('NORAD %s: animation-grid min %.4f km at t=%g s\n', ...
        noradIDs(k),min(distKm(k,:)),t(find(distKm(k,:)==min(distKm(k,:)),1)));
end
% Build the native satellite scenario using sampled C63 ephemeris.
% Ephemeris positions and velocities are in inertial coordinates (metres, m/s).
scView=satelliteScenario(startUTC,startUTC+seconds(durationS),frameStep);
satV=zeros(3,N);
for j=1:N
    if t(j)<t1
        y=deval(s1,t(j));
    elseif t(j)<t2
        y=deval(s2,t(j));
    else
        y=deval(s3,t(j));
    end
    satV(:,j)=y(4:6);
end
% A velocity discontinuity at a burn is represented at the sample instant.
missionV3=satellite(scView,satP,satV,'CoordinateFrame','inertial', ...
    'Name','C63 maneuvering satellite');
debrisV3=cell(1,numel(noradIDs));
for k=1:numel(noradIDs)
    debrisV3{k}=satellite(scView,char(tmpFiles(k)), ...
        'OrbitPropagator','sgp4','Name',char("NORAD "+noradIDs(k)));
end
viewerHandle=satelliteScenarioViewer(scView);
show(missionV3);
for k=1:numel(noradIDs),show(debrisV3{k});end
fprintf('V3 viewer ready: 10-second display samples; C63 burns preserved.\n');
fprintf('Close visual markers do not imply collision. Check numeric distances.\n');
fprintf('Starting satelliteScenario playback...\n');
play(scView);
fprintf('Done. Visualization of two selected debris only.\n');
function dy=rhs(y,mu)
r=y(1:3);dy=[y(4:6);-mu*r/norm(r)^3];
end
function T=tangent(y)
r=y(1:3);v=y(4:6);R=r/norm(r);N=cross(r,v);N=N/norm(N);
T=cross(N,R);T=T/norm(T);
end
function p=positionsAt(sols,t,t1,t2)
p=zeros(3,numel(t));i1=t<t1;i2=t>=t1&t<t2;i3=t>=t2;
if any(i1),p(:,i1)=deval(sols{1},t(i1),1:3);end
if any(i2),p(:,i2)=deval(sols{2},t(i2),1:3);end
if any(i3),p(:,i3)=deval(sols{3},t(i3),1:3);end
end
function deleteTemps(paths)
for i=1:numel(paths)
if strlength(paths(i))>0 && isfile(paths(i)),delete(paths(i));end
end
end
