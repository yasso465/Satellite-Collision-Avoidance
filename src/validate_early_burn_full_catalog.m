%% FULL CATALOG NOMINAL VALIDATION - three early-burn candidates
% MATLAB R2026b / Aerospace Toolbox
% 16,883 original 3LE objects (or count in debris_matlab.tle).
% 1-second global screening; 0.1-second local refinement for <30 km.
% Propagation errors reported, NEVER silently excluded.
% NOTE: Numerical screening only, not a continuous-time proof.
clear; clc;
tleFile='debris_matlab.tle';
assert(isfile(tleFile),'Missing debris_matlab.tle');
startTime=datetime(2026,10,12,0,0,0,'TimeZone','UTC');
missionLength=10800; mu=3.986004418e14; Re=6378137;
opts=odeset('RelTol',1e-10,'AbsTol',1e-5);
% Candidate 3, 63, 123 (same times, different second burn)
labels=["C3_050";"C63_060";"C123_070"];
burns=[2600 8000 -0.40 +0.10; ...
       2600 8000 -0.40 +0.20; ...
       2600 8000 -0.40 +0.30];
thresholdKm=10.5; refineIfBelowKm=30;
fineStep=0.1; fineHalfWindow=12;
coarseStep=1; progressEvery=100;
outCSV='early_burn_full_catalog_results.csv';
summaryCSV='early_burn_full_catalog_summary.csv';
lines=readlines(tleFile);lines=lines(strlength(strtrim(lines))>0);
assert(mod(numel(lines),3)==0,'Expected 3-line TLE');
nObj=numel(lines)/3; nC=size(burns,1);
sc=satelliteScenario(startTime,startTime+seconds(missionLength),1);
mission=satellite(sc,Re+500e3,0,51.6,0,0,0, ...
    'OrbitPropagator','two-body-keplerian');
[p0,v0]=states(mission,'CoordinateFrame','inertial');
y0=[p0(:,1);v0(:,1)]; refEl=orbitalElements(y0,mu);
tt=0:coarseStep:missionLength;
satP=cell(nC,1); sols=cell(nC,1);
for k=1:nC
    v=burns(k,:);
    [sol,yf]=propagateBurns(y0,v(1),v(2),v(3),v(4), ...
        [0 0],[0 0],mu,opts,missionLength);
    sols{k}=sol;
    satP{k}=positionsAt(sol,tt,v(1),v(2));
    el=orbitalElements(yf,mu);
    fprintf('%s: totalDV %.3f m/s | da %.6f km | de %.9f | orbitPass %d\n', ...
        labels(k),abs(v(3))+abs(v(4)),abs(el(1)-refEl(1))/1000, ...
        abs(el(2)-refEl(2)), ...
        abs(el(1)-refEl(1))<=5000 && abs(el(2)-refEl(2))<=1e-3);
end
% Preallocate per-object records; one row per object x candidate.
nRows=nObj*nC;
objectIndex=repelem((1:nObj)',nC);
candidate=repmat(labels,nObj,1);
norad=strings(nRows,1); name=strings(nRows,1);
status=repmat("NOT_RUN",nRows,1);
min1=nan(nRows,1);tca1=nan(nRows,1);
minFine=nan(nRows,1);tcaFine=nan(nRows,1);
refined=false(nRows,1);errorText=strings(nRows,1);
tic;
fprintf('\nFULL CATALOG: %d objects x %d candidates, 1 s global...\n',nObj,nC);
for i=1:nObj
    idx=(i-1)*nC+(1:nC);
    name(idx)=strtrim(lines(3*i-2));
    l1=char(lines(3*i-1));
    if numel(l1)>=7, norad(idx)=string(strtrim(l1(3:7)));end
    tmp=[tempname,'.tle']; fid=fopen(tmp,'w');
    if fid<0
        status(idx)="ERROR";errorText(idx)="Cannot create temporary TLE";
        continue
    end
    for q=1:3
        fprintf(fid,'%s\n',char(lines(3*i-3+q)));
    end
    fclose(fid);
    try
        sd=satelliteScenario(startTime,startTime+seconds(missionLength),coarseStep);
        debris=satellite(sd,tmp,'OrbitPropagator','sgp4');
        [pD,~,timesD]=states(debris,'CoordinateFrame','inertial');
        tD=seconds(timesD(:)-startTime)';
        assert(size(pD,1)==3 && size(pD,2)==numel(tt) && ...
            all(isfinite(pD(:))) && max(abs(tD-tt))<1e-6, ...
            'Global debris propagation grid mismatch or nonfinite');
        for k=1:nC
            row=idx(k);
            d=vecnorm(pD-satP{k},2,1)/1000;
            [min1(row),ix]=min(d);
            tca1(row)=tt(ix);
            status(row)="OK";
            % Fine sampling only near potentially relevant encounters.
            % Search around EVERY local 1s minimum below 30 km,
            % plus global-minimum index. Merge overlapping intervals.
            if min1(row)<refineIfBelowKm
                seeds=find(d<refineIfBelowKm);
                mask=false(size(tt));
                for seed=seeds(:)'
                    mask(max(1,seed-ceil(fineHalfWindow)): ...
                         min(numel(tt),seed+ceil(fineHalfWindow)))=true;
                end
                edges=diff([false mask false]);
                starts=find(edges==1); ends=find(edges==-1)-1;
                best=inf; bestT=nan;
                for w=1:numel(starts)
                    ta=tt(starts(w));tb=tt(ends(w));
                    % If entire mission is close, cap per-window work by
                    % splitting would be necessary; normally not needed.
                    sf=satelliteScenario(startTime+seconds(ta), ...
                        startTime+seconds(tb),fineStep);
                    df=satellite(sf,tmp,'OrbitPropagator','sgp4');
                    [pf,~,tf]=states(df,'CoordinateFrame','inertial');
                    tFine=seconds(tf(:)-startTime)';
                    ps=positionsAt(sols{k},tFine,burns(k,1),burns(k,2));
                    [md,ixf]=min(vecnorm(pf-ps,2,1)/1000);
                    if md<best,best=md;bestT=tFine(ixf);end
                end
                minFine(row)=best;tcaFine(row)=bestT;
                refined(row)=true;
            end
        end
    catch ME
        status(idx)="ERROR";
        errorText(idx)=string(ME.message);
        min1(idx)=nan;tca1(idx)=nan;
        minFine(idx)=nan;tcaFine(idx)=nan;refined(idx)=false;
    end
    if isfile(tmp),delete(tmp);end
    if mod(i,progressEvery)==0 || i==nObj
        fprintf('%d/%d | %.1f min | failures %d | refined objects %d\n', ...
            i,nObj,toc/60,nnz(status(1:i*nC)=="ERROR")/nC, ...
            nnz(refined(1:i*nC)));
        % Save progress; rerunning starts from beginning, not resume.
        if mod(i,1000)==0
            save('early_burn_full_catalog_checkpoint.mat', ...
                'i','status','min1','tca1','minFine','tcaFine', ...
                'refined','errorText','-v7.3');
        end
    end
end
R=table(objectIndex,candidate,norad,name,status,min1,tca1, ...
    minFine,tcaFine,refined,errorText, ...
    'VariableNames',{'ObjectIndex','Candidate','NORAD','Name','Status', ...
    'Min_1s_km','TCA_1s_s','Min_01s_km','TCA_01s_s', ...
    'Refined','Error'});
writetable(R,outCSV);
S=table('Size',[nC 9], ...
    'VariableTypes',{'string','double','double','double','double', ...
    'double','double','double','double'}, ...
    'VariableNames',{'Candidate','TotalDV_mps','ObjectsOK', ...
    'ObjectsError','GlobalMin_1s_km','GlobalMin_01s_km', ...
    'Below10p5_1s','Below10p5_01s','RefinedCount'});
for k=1:nC
    x=R(R.Candidate==labels(k),:);
    S.Candidate(k)=labels(k);
    S.TotalDV_mps(k)=sum(abs(burns(k,3:4)));
    S.ObjectsOK(k)=nnz(x.Status=="OK");
    S.ObjectsError(k)=nnz(x.Status=="ERROR");
    S.GlobalMin_1s_km(k)=min(x.Min_1s_km,[],'omitnan');
    S.GlobalMin_01s_km(k)=min(x.Min_01s_km,[],'omitnan');
    S.Below10p5_1s(k)=nnz(x.Min_1s_km<thresholdKm);
    S.Below10p5_01s(k)=nnz(x.Min_01s_km<thresholdKm);
    S.RefinedCount(k)=nnz(x.Refined);
end
writetable(S,summaryCSV);
fprintf('\n===== FULL CATALOG EARLY-BURN VALIDATION =====\n');disp(S);
fprintf('Ten closest objects per candidate (refined if available):\n');
for k=1:nC
    x=R(R.Candidate==labels(k) & R.Status=="OK",:);
    x.EffectiveMin_km=x.Min_1s_km;
    ix=isfinite(x.Min_01s_km);
    x.EffectiveMin_km(ix)=x.Min_01s_km(ix);
    x=sortrows(x,'EffectiveMin_km');
    fprintf('\n%s:\n',labels(k));disp(x(1:min(10,height(x)), ...
        {'ObjectIndex','NORAD','EffectiveMin_km','TCA_1s_s','TCA_01s_s'}));
end
fprintf('Outputs: %s, %s\n',outCSV,summaryCSV);
fprintf(['CAUTION: 1s global numerical sampling, 0.1s refinement for <30km, ' ...
    'not a rigorous continuous-time bound. SGP4 failures remain unresolved.\n']);

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
