%% SEARCH BURN TIMES AND MAGNITUDES - preliminary local screening
% MATLAB R2026b, Aerospace Toolbox. Run in Satellite-Collision-Avoidance.
% Searches earlier t1 to affect debris 9708 (~4624s).
% Wide local windows: +/- 300 s, 1 s SGP4. NOT full catalog validation.
% Stage 1 nominal screening; Stage 2 paired 250-trial Monte Carlo shortlist.
clear; clc;
tleFile='debris_matlab.tle'; inputFile='robust09_high_precision_validation.csv';
assert(isfile(tleFile)&&isfile(inputFile),'Missing input files');
T=readtable(inputFile,'TextType','string');
assert(all(ismember({'ObjectIndex','Status','RefinedTCA_s','RefinedMin_km'}, ...
    T.Properties.VariableNames)),'Unexpected CSV columns');
assert(all(T.Status=="OK"),'Failed local encounters');
T=sortrows(T,'RefinedMin_km');
startTime=datetime(2026,10,12,0,0,0,'TimeZone','UTC');
missionLength=10800; mu=3.986004418e14; Re=6378137;
opts=odeset('RelTol',1e-10,'AbsTol',1e-5);
sc=satelliteScenario(startTime,startTime+seconds(missionLength),1);
mission=satellite(sc,Re+500e3,0,51.6,0,0,0, ...
    'OrbitPropagator','two-body-keplerian');
[p0,v0]=states(mission,'CoordinateFrame','inertial');
y0=[p0(:,1);v0(:,1)]; refEl=orbitalElements(y0,mu);
lines=readlines(tleFile);lines=lines(strlength(strtrim(lines))>0);
assert(mod(numel(lines),3)==0,'Expected 3-line TLE');
assert(all(T.ObjectIndex>=1 & T.ObjectIndex<=numel(lines)/3),'Bad object index');
% Grid: 5 * 3 * 4 * 3 = 180 + baseline
t1Grid=[2000 2600 3200 3800 4400];
t2Grid=[8000 8800 9600];
dv1Grid=[-0.4 -0.6 -0.8 -1.0];
dv2Grid=[0.1 0.2 0.3];
[a,b,c,d]=ndgrid(t1Grid,t2Grid,dv1Grid,dv2Grid);
C=unique([5489 8789 -.70 .20; a(:) b(:) c(:) d(:)],'rows','stable');
nC=size(C,1); nEnc=height(T); threshold=10.5;
halfWindow=300; step=1;
E=struct('id',{},'tt',{},'pD',{},'tca',{});
fprintf('Caching %d SGP4 objects, +/- %d s at %.1f s...\n', ...
    nEnc,halfWindow,step);
for i=1:nEnc
    k=T.ObjectIndex(i); tc=double(T.RefinedTCA_s(i));
    ta=max(0,floor(tc-halfWindow)); tb=min(missionLength,ceil(tc+halfWindow));
    sd=satelliteScenario(startTime+seconds(ta),startTime+seconds(tb),step);
    f=[tempname '.tle']; fid=fopen(f,'w');assert(fid>=0,'Cannot write temp TLE');
    for q=1:3, fprintf(fid,'%s\n',char(lines(3*k-3+q)));end
    fclose(fid);
    try
        sat=satellite(sd,f,'OrbitPropagator','sgp4');
        [pd,~,td]=states(sat,'CoordinateFrame','inertial');
    catch ME
        if isfile(f),delete(f);end
        rethrow(ME);
    end
    if isfile(f),delete(f);end
    tt=seconds(td(:)-startTime)';
    assert(all(isfinite(pd(:))) && size(pd,2)==numel(tt),'SGP4 cache invalid');
    E(i).id=k;E(i).tt=tt;E(i).pD=pd;E(i).tca=tc;
end
% Stage 1: nominal, no uncertainty.
M=nan(nC,10);
timer=tic;
fprintf('Stage 1: nominal screening of %d candidate burns...\n',nC);
for k=1:nC
    v=C(k,:);
    [sol,yf]=propagateBurns(y0,v(1),v(2),v(3),v(4), ...
        [0 0],[0 0],mu,opts,missionLength);
    md=inf;worst=nan;tWorst=nan;atBoundary=false;
    for i=1:nEnc
        e=E(i); ps=positionsAt(sol,e.tt,v(1),v(2));
        [di,ix]=min(vecnorm(e.pD-ps,2,1)/1000);
        if di<md
            md=di;worst=e.id;tWorst=e.tt(ix);
            atBoundary=(ix==1 || ix==numel(e.tt)) && ...
                e.tt(ix)~=0 && e.tt(ix)~=missionLength;
        end
    end
    el=orbitalElements(yf,mu);
    da=abs(el(1)-refEl(1))/1000; de=abs(el(2)-refEl(2));
    M(k,:)=[k,sum(abs(v(3:4))),md,worst,tWorst,da,de, ...
        double(da<=5 && de<=1e-3),double(atBoundary),0];
    if mod(k,20)==0 || k==nC
        fprintf('%d/%d | %.1f min\n',k,nC,toc(timer)/60);
    end
end
N=array2table([C M(:,1:9)],'VariableNames', ...
    {'T1_s','T2_s','DV1_mps','DV2_mps','Candidate','TotalDV_mps', ...
    'NominalMin_1s_km','WorstObject','TCA_s','DeltaA_km','DeltaE', ...
    'OrbitPass','WorstAtWindowBoundary'});
writetable(N,'early_burn_nominal_all.csv');
% Select viable solutions: nominal 1-s >= 10.8 km, orbital pass, no boundary.
viable=N(N.OrbitPass==1 & N.NominalMin_1s_km>=10.8 & ...
    N.WorstAtWindowBoundary==0,:);
viable=sortrows(viable,{'NominalMin_1s_km','TotalDV_mps'}, ...
    {'descend','ascend'});
nShort=min(8,height(viable));
fprintf('Nominal viable candidates: %d; shortlist: %d\n',height(viable),nShort);
if nShort==0
    fprintf('No viable candidate. Expand grid or relax screening, do NOT claim safety.\n');
    return
end
% Include reference B for paired comparison (even if outside viable).
S=viable(1:nShort,:);
baseline=N(N.Candidate==1,:);
if ~ismember(1,S.Candidate),S=[baseline;S];end
% Stage 2: exactly paired perturbations for each candidate.
nTrials=250; rng(20261016,'twister');
zT1=randn(nTrials,1);zT2=randn(nTrials,1);
zD1=randn(nTrials,1);zD2=randn(nTrials,1);
zC1=randn(nTrials,2);zC2=randn(nTrials,2);
zP=randn(3,nEnc,nTrials);zV=randn(3,nEnc,nTrials);
sigTime=1;sigDV=.01;sigCross=.005;sigP=250;sigV=.05;
Q=nan(height(S),9);
fprintf('Stage 2: %d candidates x %d paired trials...\n',height(S),nTrials);
for k=1:height(S)
    v=[S.T1_s(k),S.T2_s(k),S.DV1_mps(k),S.DV2_mps(k)];
    mins=nan(nTrials,1);passes=false(nTrials,1);
    for j=1:nTrials
        t1=v(1)+sigTime*zT1(j);t2=v(2)+sigTime*zT2(j);
        [sol,yf]=propagateBurns(y0,t1,t2, ...
            v(3)+sigDV*zD1(j),v(4)+sigDV*zD2(j), ...
            sigCross*zC1(j,:),sigCross*zC2(j,:), ...
            mu,opts,missionLength);
        md=inf;
        for i=1:nEnc
            e=E(i); tt=e.tt;
            pD=e.pD+sigP*zP(:,i,j)+sigV*zV(:,i,j)*(tt-e.tca);
            ps=positionsAt(sol,tt,t1,t2);
            md=min(md,min(vecnorm(pD-ps,2,1)/1000));
        end
        mins(j)=md;
        el=orbitalElements(yf,mu);
        passes(j)=abs(el(1)-refEl(1))<=5000 && ...
            abs(el(2)-refEl(2))<=1e-3;
    end
    q=prctile(mins,[1 5]);
    Q(k,:)=[S.Candidate(k),nnz(mins<threshold),min(mins), ...
        q,nnz(passes),S.TotalDV_mps(k),v(1),v(2)];
    fprintf('Candidate %d | t1=%g t2=%g | dv=[%+.2f %+.2f] | breaches %d/%d | min %.4f km | orbit %d/%d\n', ...
        S.Candidate(k),v,nnz(mins<threshold),nTrials,min(mins),nnz(passes),nTrials);
end
R=array2table(Q,'VariableNames', ...
    {'Candidate','Breaches_1s','Minimum_1s_km','P01_1s_km', ...
    'P05_1s_km','OrbitPassCount','TotalDV_mps','T1_s','T2_s'});
R=sortrows(R,{'OrbitPassCount','Breaches_1s','P01_1s_km', ...
    'TotalDV_mps'},{'descend','ascend','descend','ascend'});
writetable(R,'early_burn_monte_carlo_shortlist.csv');
fprintf('\n===== EARLY BURN TIME SEARCH =====\n');disp(R);
fprintf('Files: early_burn_nominal_all.csv, early_burn_monte_carlo_shortlist.csv\n');
fprintf(['WARNING: only 12 historical encounter windows, 1-s grid. ' ...
    'New burns can shift encounters OUTSIDE these windows or create NEW encounters. ' ...
    'All candidates require full-catalog screening and 0.1-s refinement.\n']);

function [sol,yf]=propagateBurns(y0,t1,t2,dv1,dv2,cross1,cross2,mu,opts,tEnd)
    s1=ode113(@(t,y)twoBody(t,y,mu),[0 t1],y0,opts);
    y1=deval(s1,t1);y1(4:6)=y1(4:6)+burnVector(y1,dv1,cross1);
    s2=ode113(@(t,y)twoBody(t,y,mu),[t1 t2],y1,opts);
    y2=deval(s2,t2);y2(4:6)=y2(4:6)+burnVector(y2,dv2,cross2);
    s3=ode113(@(t,y)twoBody(t,y,mu),[t2 tEnd],y2,opts);
    yf=deval(s3,tEnd);sol={s1,s2,s3};
end
function d=burnVector(y,tangentialDV,crossDV)
    r=y(1:3);v=y(4:6);R=r/norm(r);
    N=cross(r,v);N=N/norm(N);T=cross(N,R);T=T/norm(T);
    d=tangentialDV*T+crossDV(1)*R+crossDV(2)*N;
end
function pos=positionsAt(sol,t,t1,t2)
    pos=zeros(3,numel(t));
    i1=t<t1;i2=t>=t1 & t<t2;i3=t>=t2;
    if any(i1),pos(:,i1)=deval(sol{1},t(i1),1:3);end
    if any(i2),pos(:,i2)=deval(sol{2},t(i2),1:3);end
    if any(i3),pos(:,i3)=deval(sol{3},t(i3),1:3);end
end
function dy=twoBody(~,y,mu)
    r=y(1:3);dy=[y(4:6);-mu*r/norm(r)^3];
end
function el=orbitalElements(y,mu)
    r=y(1:3);v=y(4:6);rn=norm(r);vn=norm(v);
    a=1/(2/rn-vn^2/mu);
    ev=((vn^2-mu/rn)*r-dot(r,v)*v)/mu;
    el=[a;norm(ev)];
end

function d=localDistance(t,tt,pD,dP,dV,tRef,sol,t1,t2)
    % pchip interpolates each Cartesian component of the cached SGP4 orbit.
    t=t(:)';
    pNom=interp1(tt,pD.',t,'pchip').';
    pPert=pNom+dP+dV*(t-tRef);
    pSat=positionsAt(sol,t,t1,t2);
    d=vecnorm(pPert-pSat,2,1)/1000;
    if isscalar(d),d=double(d);end
end
