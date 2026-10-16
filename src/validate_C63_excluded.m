%% C63: validation of excluded TLEs against maneuvered satellite
% Requires debris_leo.tle (original) and debris_filtered.tle (original IDs).
% Resume-safe CSV. MATLAB R2026b + Aerospace Toolbox.
clear;clc;
originalFile='debris_leo.tle'; filteredFile='debris_filtered.tle';
outputFile='C63_excluded_validation.csv';
assert(isfile(originalFile)&&isfile(filteredFile), ...
    'Place debris_leo.tle et debris_filtered.tle dans le dossier courant.');
orig=readlines(originalFile);orig=orig(strlength(strtrim(orig))>0);
filt=readlines(filteredFile);filt=filt(strlength(strtrim(filt))>0);
assert(mod(numel(orig),3)==0 && mod(numel(filt),3)==0,'3LE invalid');
nOrig=numel(orig)/3;nFilt=numel(filt)/3;
okey=orig(2:3:end)+"|"+orig(3:3:end);
fkey=filt(2:3:end)+"|"+filt(3:3:end);
assert(numel(unique(okey))==nOrig && numel(unique(fkey))==nFilt, ...
    'Duplicate TLE entries: ambiguous catalog mapping');
assert(all(ismember(fkey,okey)),'Filtered catalog does not match original');
excluded=find(~ismember(okey,fkey));
assert(numel(excluded)==nOrig-nFilt,'Catalog counts inconsistent');
fprintf('Original %d | retained %d | excluded %d\n',nOrig,nFilt,numel(excluded));
ids=string(extractBetween(orig(2:3:end),3,7));
numMask=~cellfun(@isempty,regexp(cellstr(ids),'^[0-9]{5}$','once'));
used=unique(str2double(ids(numMask)));
freeIds=setdiff((0:99999)',used);
alpha=excluded(~numMask(excluded));
assert(numel(freeIds)>=numel(alpha),'No free temporary IDs');
replacement=strings(nOrig,1);
for k=1:numel(alpha)
    replacement(alpha(k))=string(sprintf('%05d',freeIds(k)));
end
startTime=datetime(2026,10,12,0,0,0,'TimeZone','UTC');
tEnd=10800;mu=3.986004418e14;Re=6378137;
t1=2600;t2=8000;dv1=-0.40;dv2=+0.20;
threshold=10.5;refineBelow=40;fineStep=0.1;halfWindow=12;
opts=odeset('RelTol',1e-10,'AbsTol',1e-5);
sc=satelliteScenario(startTime,startTime+seconds(tEnd),1);
sat=satellite(sc,Re+500e3,0,51.6,0,0,0,'OrbitPropagator','two-body-keplerian');
[p0,v0]=states(sat,'CoordinateFrame','inertial');y0=[p0(:,1);v0(:,1)];
[sol,~]=propagateBurns(y0,t1,t2,dv1,dv2,[0 0],[0 0],mu,opts,tEnd);
tt=0:tEnd;satP=positionsAt(sol,tt,t1,t2);
if isfile(outputFile)
    old=readtable(outputFile,'TextType','string');
    assert(all(ismember(['OriginalIndex','Status'],string(old.Properties.VariableNames))), ...
        'Existing output has wrong schema');
    assert(all(ismember(old.OriginalIndex,excluded)), ...
        'Existing output belongs to a different catalog');
    done=unique(old.OriginalIndex);
else
    done=[];
end
pending=setdiff(excluded,done,'stable');
fprintf('Previously saved %d | pending %d\n',numel(done),numel(pending));
rows=table('Size',[0 9], ...
    'VariableTypes',{'double','string','string','double','double', ...
                     'double','double','logical','string'}, ...
    'VariableNames',{'OriginalIndex','NORAD','Status','Min_1s_km', ...
        'TCA_1s_s','Min_01s_km','TCA_01s_s','Refined','Error'});
tic;
for j=1:numel(pending)
    i=pending(j); tempFile=[tempname,'.tle'];
    min1=nan;tca1=nan;minFine=nan;tcaFine=nan;
    refined=false;status="OK";err="";
    try
        triplet=orig(3*i-2:3*i);
        if replacement(i)~=""
            for q=2:3
                line=char(triplet(q));
                assert(numel(line)>=69,'Invalid Alpha-5 TLE length');
                line(3:7)=char(replacement(i));
                line(69)=tleChecksumDigit(line);
                triplet(q)=string(line);
            end
        end
        fid=fopen(tempFile,'w');assert(fid>=0,'Cannot create temporary TLE');
        for q=1:3,fprintf(fid,'%s\n',char(triplet(q)));end
        fclose(fid);
        sd=satelliteScenario(startTime,startTime+seconds(tEnd),1);
        debris=satellite(sd,tempFile,'OrbitPropagator','sgp4');
        [pd,~,times]=states(debris,'CoordinateFrame','inertial');
        timeS=seconds(times(:)-startTime)';
        assert(size(pd,1)==3 && size(pd,2)==numel(tt) && ...
            all(isfinite(pd(:))) && max(abs(timeS-tt))<1e-6, ...
            'SGP4 output invalid or time mismatch');
        d=vecnorm(pd-satP,2,1)/1000;
        [min1,ix]=min(d);tca1=tt(ix);
        if min1<refineBelow
            seeds=find(d<refineBelow);
            mask=false(size(tt));
            for seed=seeds(:)'
                mask(max(1,seed-halfWindow):min(numel(tt),seed+halfWindow))=true;
            end
            edges=diff([false mask false]);starts=find(edges==1);
            ends=find(edges==-1)-1;
            for w=1:numel(starts)
                ta=tt(starts(w));tb=tt(ends(w));
                sf=satelliteScenario(startTime+seconds(ta), ...
                    startTime+seconds(tb),fineStep);
                df=satellite(sf,tempFile,'OrbitPropagator','sgp4');
                [pf,~,tf]=states(df,'CoordinateFrame','inertial');
                tfine=seconds(tf(:)-startTime)';
                ps=positionsAt(sol,tfine,t1,t2);
                [md,ii]=min(vecnorm(pf-ps,2,1)/1000);
                if isnan(minFine)||md<minFine
                    minFine=md;tcaFine=tfine(ii);
                end
            end
            refined=true;
        end
    catch ME
        status="ERROR";err=string(ME.message);
        fprintf('SGP4 ERROR original index %d, NORAD %s: %s\n', ...
            i,ids(i),ME.message);
    end
    if isfile(tempFile),delete(tempFile);end
    rows(end+1,:)={i,ids(i),status,min1,tca1,minFine,tcaFine,refined,err}; %#ok<AGROW>
    if mod(j,25)==0 || j==numel(pending)
        if isfile(outputFile)
            writetable(rows,outputFile,'WriteMode','append');
        else
            writetable(rows,outputFile);
        end
        rows(1:height(rows),:)=[];
    end
    if mod(j,100)==0 || j==numel(pending)
        fprintf('%d/%d | %.1f min\n',j,numel(pending),toc/60);
    end
end
R=readtable(outputFile,'TextType','string');
[~,ia]=unique(flipud(R.OriginalIndex),'stable');
R=R(sort(height(R)-ia+1),:);
assert(all(ismember(R.OriginalIndex,excluded)),'Unexpected indices');
valid=R.Status=="OK";
effective=R.Min_1s_km;
fine=isfinite(R.Min_01s_km);effective(fine)=R.Min_01s_km(fine);
R.EffectiveMin_km=effective;
R=sortrows(R,'EffectiveMin_km');
fprintf('\n===== C63 EXCLUDED CATALOG VALIDATION =====\n');
fprintf('Expected excluded: %d | completed: %d\n',numel(excluded),height(R));
fprintf('SGP4 OK: %d | failures: %d\n',nnz(valid),nnz(~valid));
fprintf('Below 10.5 km at 1s: %d | after refinement: %d\n', ...
    nnz(valid & R.Min_1s_km<threshold),nnz(valid & effective<threshold));
fprintf('Smallest observed separation: %.6f km\n',min(effective,[],'omitnan'));
fprintf('Closest 10 objects:\n');
disp(R(1:min(10,height(R)),{'OriginalIndex','NORAD','Status', ...
    'EffectiveMin_km','TCA_1s_s','TCA_01s_s'}));
if any(~valid)
    fprintf('SGP4 failure details:\n');
    disp(R(~valid,{'OriginalIndex','NORAD','Error'}));
end
fprintf('CSV: %s\n',outputFile);
fprintf('WARNING: finite-step numerical screening, not a continuous-time guarantee.\n');

function digit=tleChecksumDigit(line)
    chars=line(1:68);v=double(chars)-double('0');
    mask=chars>='0' & chars<='9';
    digit=char(double('0')+mod(sum(v(mask))+sum(chars=='-'),10));
end

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
