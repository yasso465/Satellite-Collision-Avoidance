%% PAIRED 10,000-TRIAL MONTE CARLO: C3 (0.50) vs C63 (0.60)
% MATLAB R2026b / Aerospace Toolbox.
% Requires debris_matlab.tle and early_burn_full_catalog_results.csv
% Run from project directory. Local risk screening, NOT collision probability.
clear;clc;
fileTLE='debris_matlab.tle';fileCSV='early_burn_full_catalog_results.csv';
assert(isfile(fileTLE)&&isfile(fileCSV),'Missing TLE or full-catalog CSV');
N=10000; rng(20261016,'twister'); thresholdKm=10.5;
sigmaP=250;sigmaV=0.05;sigmaTime=1;sigmaTang=0.01;sigmaCross=0.005;
windowHalf=45;cacheStep=0.1;screenStep=1; refineMarginKm=2;
missionLength=10800;mu=3.986004418e14;Re=6378137;
startTime=datetime(2026,10,12,0,0,0,'TimeZone','UTC');
opts=odeset('RelTol',1e-10,'AbsTol',1e-5);
sc=satelliteScenario(startTime,startTime+seconds(missionLength),1);
mission=satellite(sc,Re+500e3,0,51.6,0,0,0,'OrbitPropagator','two-body-keplerian');
[p0,v0]=states(mission,'CoordinateFrame','inertial');
y0=[p0(:,1);v0(:,1)];refEl=orbitalElements(y0,mu);
lines=readlines(fileTLE);lines=lines(strlength(strtrim(lines))>0);
assert(mod(numel(lines),3)==0,'Expected 3LE TLE');
T=readtable(fileCSV,'TextType','string');
needed={'ObjectIndex','Candidate','Status','Min_1s_km','TCA_1s_s','Min_01s_km','TCA_01s_s'};
assert(all(ismember(needed,T.Properties.VariableNames)),'Unexpected CSV columns');
assert(all(T.Status=="OK"),'Full-catalog CSV contains propagation failures');
assert(all(ismember(["C3_050","C63_060"],unique(T.Candidate))),'Missing candidates');
% Include every object within 30 km in either full-catalog nominal run.
critical=T(T.Min_1s_km<30 & ismember(T.Candidate,["C3_050","C63_060"]),:);
ids=unique(critical.ObjectIndex,'stable');
assert(~isempty(ids),'No critical objects found');
assert(all(ids>=1 & ids<=numel(lines)/3),'Object index outside TLE');
% For each object, use union of C3/C63 nominal TCA windows.
E=struct('id',{},'tt',{},'pd',{},'tc',{});
fprintf('Caching %d critical objects (union across both candidates)...\n',numel(ids));
for ii=1:numel(ids)
    id=ids(ii); r=critical(critical.ObjectIndex==id,:);
    centers=r.TCA_01s_s;bad=~isfinite(centers);centers(bad)=r.TCA_1s_s(bad);
    for ic=1:numel(centers)
        tc=centers(ic); ta=max(0,floor(tc-windowHalf));tb=min(missionLength,ceil(tc+windowHalf));
        % Merge identical/overlapping windows by skipping if already covered.
        duplicate=false;
        for j=1:numel(E)
            if E(j).id==id && ta>=E(j).tt(1) && tb<=E(j).tt(end)
                duplicate=true;break
            end
        end
        if duplicate,continue;end
        tmp=[tempname '.tle'];fid=fopen(tmp,'w');assert(fid>=0,'Cannot create TLE');
        for q=1:3,fprintf(fid,'%s\n',char(lines(3*id-3+q)));end
        fclose(fid);
        try
            sd=satelliteScenario(startTime+seconds(ta),startTime+seconds(tb),cacheStep);
            d=satellite(sd,tmp,'OrbitPropagator','sgp4');
            [pd,~,dt]=states(d,'CoordinateFrame','inertial');
            tt=seconds(dt(:)-startTime)';
            assert(size(pd,1)==3 && size(pd,2)==numel(tt) && all(isfinite(pd(:))), ...
                'Invalid cached SGP4 states');
        catch ME
            if isfile(tmp),delete(tmp);end
            rethrow(ME);
        end
        if isfile(tmp),delete(tmp);end
        E(end+1).id=id;E(end).tt=tt;E(end).pd=pd;E(end).tc=tc; %#ok<SAGROW>
    end
    if mod(ii,10)==0 || ii==numel(ids),fprintf('Cached %d/%d objects\n',ii,numel(ids));end
end
nEnc=numel(E);fprintf('Total local windows: %d\n',nEnc);
% Paired errors for burns and each object. Same object gets same errors
% even if its nominal encounters occupy two windows.
zT=randn(N,2);zD=randn(N,2);zC=randn(N,4);
zP=randn(3,numel(ids),N);zV=randn(3,numel(ids),N);
[~,idMap]=ismember([E.id],ids);assert(all(idMap>0));
labels=["C3_050" "C63_060"];
burns=[2600 8000 -0.4 0.1;2600 8000 -0.4 0.2];
resMin=nan(N,2);resWorst=nan(N,2);resTCA=nan(N,2);
resDa=nan(N,2);resDe=nan(N,2);resPass=false(N,2);
for k=1:2
    b=burns(k,:);tic;fprintf('\nRunning %s: %d trials...\n',labels(k),N);
    for j=1:N
        t1=b(1)+sigmaTime*zT(j,1);t2=b(2)+sigmaTime*zT(j,2);
        [sol,yf]=propagateBurns(y0,t1,t2, ...
            b(3)+sigmaTang*zD(j,1),b(4)+sigmaTang*zD(j,2), ...
            sigmaCross*zC(j,1:2),sigmaCross*zC(j,3:4),mu,opts,missionLength);
        el=orbitalElements(yf,mu);
        resDa(j,k)=abs(el(1)-refEl(1))/1000;
        resDe(j,k)=abs(el(2)-refEl(2));
        resPass(j,k)=resDa(j,k)<=5 && resDe(j,k)<=1e-3;
        coarseM=inf(nEnc,1);coarseT=nan(nEnc,1);
        % First screen each local window at 1-second sampling.
        for e=1:nEnc
            ce=E(e);ix=1:round(screenStep/cacheStep):numel(ce.tt);
            tt=ce.tt(ix);p=ce.pd(:,ix);
            obj=idMap(e);
            p=p+sigmaP*zP(:,obj,j)+sigmaV*zV(:,obj,j)*(tt-ce.tc);
            ps=positionsAt(sol,tt,t1,t2);
            [coarseM(e),ii]=min(vecnorm(p-ps,2,1)/1000);
            coarseT(e)=tt(ii);
        end
        % Refine all windows close to the worst screened approach.
        [best,~]=min(coarseM);bestId=nan;bestT=nan;
        ixList=find(coarseM<=best+refineMarginKm);
        for e=ixList(:)'
            ce=E(e);tc=coarseT(e);
            lo=max(ce.tt(1),tc-2);hi=min(ce.tt(end),tc+2);
            if lo>=hi,continue;end
            obj=idMap(e);dp=sigmaP*zP(:,obj,j);dv=sigmaV*zV(:,obj,j);
            fun=@(t)localDistance(t,ce.tt,ce.pd,dp,dv,ce.tc,sol,t1,t2);
            % Search around 1s minimum; evaluate exact boundary as well.
            [tOpt,mOpt]=fminbnd(fun,lo,hi,optimset('TolX',1e-4));
            vv=[mOpt,fun(lo),fun(hi)];tv=[tOpt,lo,hi];
            [mOpt,ii]=min(vv);tOpt=tv(ii);
            if mOpt<best,best=mOpt;bestId=ce.id;bestT=tOpt;end
        end
        % If refined minima are >= screened min, identify worst screened object.
        if isnan(bestId)
            [best,ie]=min(coarseM);bestId=E(ie).id;bestT=coarseT(ie);
        end
        resMin(j,k)=best;resWorst(j,k)=bestId;resTCA(j,k)=bestT;
        if mod(j,500)==0
            fprintf('%s: %d/%d | %.1f min | breaches %d | orbit passes %d\n', ...
                labels(k),j,N,toc/60,nnz(resMin(1:j,k)<thresholdKm),nnz(resPass(1:j,k)));
            save('early_burn_mc_checkpoint.mat','j','k','resMin','resWorst', ...
                'resTCA','resDa','resDe','resPass','-v7.3');
        end
    end
end
trial=(1:N)';R=table(trial,resMin(:,1),resMin(:,2),resWorst(:,1), ...
    resWorst(:,2),resTCA(:,1),resTCA(:,2),resDa(:,1),resDa(:,2), ...
    resDe(:,1),resDe(:,2),resPass(:,1),resPass(:,2), ...
    'VariableNames',{'Trial','C3_Min_km','C63_Min_km', ...
    'C3_WorstObject','C63_WorstObject','C3_TCA_s','C63_TCA_s', ...
    'C3_DeltaA_km','C63_DeltaA_km','C3_DeltaE','C63_DeltaE', ...
    'C3_OrbitPass','C63_OrbitPass'});
writetable(R,'early_burn_mc_10000_trials.csv');
S=table('Size',[2 9],'VariableTypes', ...
    {'string','double','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'Candidate','DV_mps','Breaches10p5','Rate_percent', ...
    'Min_km','P01_km','P05_km','OrbitPassCount','Trials'});
for k=1:2
    S.Candidate(k)=labels(k);S.DV_mps(k)=sum(abs(burns(k,3:4)));
    S.Breaches10p5(k)=nnz(resMin(:,k)<thresholdKm);
    S.Rate_percent(k)=100*S.Breaches10p5(k)/N;
    S.Min_km(k)=min(resMin(:,k));qq=prctile(resMin(:,k),[1 5]);
    S.P01_km(k)=qq(1);S.P05_km(k)=qq(2);
    S.OrbitPassCount(k)=nnz(resPass(:,k));S.Trials(k)=N;
end
writetable(S,'early_burn_mc_10000_summary.csv');
fprintf('\n===== PAIRED 10000-TRIAL MONTE CARLO =====\n');disp(S);
fprintf('Paired comparison C63 > C3 separation in %.2f%% of trials\n', ...
    100*mean(resMin(:,2)>resMin(:,1)));
fprintf('Saved early_burn_mc_10000_trials.csv and early_burn_mc_10000_summary.csv\n');
fprintf(['IMPORTANT: local windows only; assumed independent Gaussian errors; ' ...
    'threshold breaches are NOT physical collision probabilities.\n']);

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
