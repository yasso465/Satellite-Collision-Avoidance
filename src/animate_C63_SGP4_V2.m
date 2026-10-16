%% C63 V2 — SGP4 scientific animation with close-encounter zoom
% MATLAB R2026b + Aerospace Toolbox. Paths resolved from script location.
% Reads actual TLEs; no synthetic debris trajectories.
% Quick preview: set exportMP4=false. Rendering may take a few minutes.
clear; clc; close all;
exportMP4 = false; fps = 24; frameStep = 5; % seconds per overview frame
% Set exportMP4=true after the preview works.
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
% Five-second overview frames are not a substitute for 0.1 s TCA refinement.
% Fine SGP4 sampling in local encounter windows provides the precise insets.
fineStep=0.1; tcaApprox=[8188.5,10796]; fineT=cell(1,numel(noradIDs)); fineD=cell(1,numel(noradIDs));
for k=1:numel(noradIDs)
    a=max(0,floor(tcaApprox(k)-15)); b=min(durationS,ceil(tcaApprox(k)+15));
    sf=satelliteScenario(startUTC+seconds(a),startUTC+seconds(b),fineStep);
    df=satellite(sf,char(tmpFiles(k)),'OrbitPropagator','sgp4');
    [pf,~,tf]=states(df,'CoordinateFrame','inertial');
    fineT{k}=seconds(tf(:)-startUTC)';
    ps=positionsAt(sols,fineT{k},t1,t2);
    fineD{k}=vecnorm(pf-ps,2,1)/1000;
    [dm,im]=min(fineD{k});
    fprintf('NORAD %s: fine-grid min %.6f km at t=%.1f s\\n',noradIDs(k),dm,fineT{k}(im));
end
% Layout: inertial overview and encounter-distance histories.
fig=figure('Color','w','Position',[80 80 1300 650]);
ax1=subplot(1,2,1,'Parent',fig);hold(ax1,'on');axis(ax1,'equal');grid(ax1,'on');
[X,Y,Z]=sphere(32);
surf(ax1,Re*X/1e3,Re*Y/1e3,Re*Z/1e3,'FaceColor',[0.14 0.36 0.66], ...
    'FaceAlpha',0.85,'EdgeColor','none');
plot3(ax1,satP(1,:)/1e3,satP(2,:)/1e3,satP(3,:)/1e3, ...
    'Color',[0.3 0.6 0.95],'LineStyle',':');
satMarker=plot3(ax1,nan,nan,nan,'o','MarkerFaceColor',[1 .65 0], ...
    'MarkerEdgeColor','k','MarkerSize',9);
colors=[.9 .2 .2;.5 .15 .85]; debMarkers=gobjects(numel(noradIDs),1);
for k=1:numel(noradIDs)
    plot3(ax1,squeeze(debP(1,:,k))/1e3,squeeze(debP(2,:,k))/1e3, ...
        squeeze(debP(3,:,k))/1e3,'Color',colors(k,:),'LineStyle','--');
    debMarkers(k)=plot3(ax1,nan,nan,nan,'o','MarkerFaceColor',colors(k,:), ...
        'MarkerEdgeColor','k','MarkerSize',7);
end
xlabel(ax1,'ECI x (km)');ylabel(ax1,'ECI y (km)');zlabel(ax1,'ECI z (km)');
view(ax1,35,24);title(ax1,'C63 orbit and SGP4 debris');
fullLimits=[-8200 8200]; xlim(ax1,fullLimits);ylim(ax1,fullLimits);zlim(ax1,fullLimits);
ax2=subplot(1,2,2,'Parent',fig);hold(ax2,'on');grid(ax2,'on');
for k=1:numel(noradIDs)
    plot(ax2,t/60,distKm(k,:),'Color',colors(k,:),'LineWidth',1.5, ...
        'DisplayName',"NORAD "+noradIDs(k));
end
yline(ax2,10.5,'k--','10.5 km screening threshold');
xline(ax2,t1/60,'--','Burn 1');xline(ax2,t2/60,'--','Burn 2');
cursor=xline(ax2,0,'k-','LineWidth',1.3,'HandleVisibility','off');
xlabel(ax2,'Mission time (min)');ylabel(ax2,'Separation (km)');
title(ax2,'Separation history (5 s overview samples)');
% Overlay independently propagated 0.1 s SGP4 local windows.
for k=1:numel(noradIDs)
    plot(ax2,fineT{k}/60,fineD{k},'Color',colors(k,:), ...
        'LineWidth',2.3,'HandleVisibility','off');
end
xlim(ax2,[0 durationS/60]);
legend(ax2,'Location','best');set(ax2,'YScale','log');
video=[];
if exportMP4
    animationDir=fullfile(projectRoot,'animations');
    if ~isfolder(animationDir),mkdir(animationDir);end
    video=VideoWriter(fullfile(animationDir,'C63_SGP4_encounters_V2.mp4'),'MPEG-4');
    video.FrameRate=fps;open(video);
end
for j=1:N
    set(satMarker,'XData',satP(1,j)/1e3,'YData',satP(2,j)/1e3, ...
        'ZData',satP(3,j)/1e3);
    for k=1:numel(noradIDs)
        set(debMarkers(k),'XData',debP(1,j,k)/1e3, ...
            'YData',debP(2,j,k)/1e3,'ZData',debP(3,j,k)/1e3);
    end
    cursor.Value=t(j)/60;
    % Automatic scientific close-up near the two actual encounter windows.
    nearIdx=find(abs(t(j)-tcaApprox)<65,1);
    if ~isempty(nearIdx)
        k=nearIdx;
        midpoint=(satP(:,j)+debP(:,j,k))/2000; % ECI km
        span=450; xlim(ax1,midpoint(1)+[-span span]);
        ylim(ax1,midpoint(2)+[-span span]);
        zlim(ax1,midpoint(3)+[-span span]);
        xlim(ax2,(tcaApprox(k)+[-80 80])/60);
        ylim(ax2,[8 2000]);
        title(ax1,sprintf('CLOSE APPROACH | NORAD %s | 0.1 s refinement',noradIDs(k)));
    else
        xlim(ax1,fullLimits);ylim(ax1,fullLimits);zlim(ax1,fullLimits);
        xlim(ax2,[0 durationS/60]);ylim(ax2,'auto');
        title(ax1,'C63 orbit and SGP4 debris');
    end
    if t(j)<t1,phase='Before burn 1';
    elseif t(j)<t2,phase='Between burns';
    else,phase='After burn 2';end
    sgtitle(fig,sprintf('C63 | %02d:%02d:%02d UTC | %s', ...
        floor(t(j)/3600),mod(floor(t(j)/60),60),mod(t(j),60),phase));
    drawnow;
    if exportMP4,writeVideo(video,getframe(fig));end
end
if exportMP4,close(video);fprintf('Saved to %s\\n',fullfile(projectRoot,'animations','C63_SGP4_encounters_V2.mp4'));end
fprintf('Done. Selected-object visualization only; not full-catalog revalidation.\n');
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
