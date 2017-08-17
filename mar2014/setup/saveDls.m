function [Dwd, Dwdt, tp, yil, dtlims, core_order, bv, Pm, d18mv, d18msv,LGMmv] = mkDls

clear

path1 = '../../data/benthic/';
lims = [-25000,-5000];
[D,tD,Di,core_order] = mkDni(path1,[-26000,-4000]); % larger limits before interp
M = max(Di);
TSTEP = 200;

D1 = [];
Pm = [];
d18mv = [];
d18msv = [];
dtlims = [];

for ii = 1:M
    ly = sum(Di==ii);
    y = D(Di==ii);
    t = tD(Di==ii);
    N = 0.2^2*eye(ly);
    
    t0 = ceil(min(t/TSTEP))*TSTEP;
    tp = t0:TSTEP:max(t);
    
    tp(tp<lims(1) | tp>lims(2)) = [];
    
    % interpolate via lsinterp
    [d18int,P,S,d18m,d18ms] = lsinterp(y,t,tp,N);
    
    d18mv = [d18mv;d18m];
    d18msv = [d18msv;d18ms];
    
    if any(P(:)<0)
        disp('Warning: some values of P<0')
    end
    P(P(:)<0) = 0;
    
    % Concatenate this record onto the D vector
    D1 = [D1;d18int(:)];
    
    % Concatenate the uncertainty matrix
    [lpmr,lpmc] = size(Pm);
    [lmr,lmc] = size(P); % square
    Pm = [[Pm,zeros(lpmr,lmc)];[zeros(lmr,lpmc),P]];
    
    dtlim = [min(tp),max(tp)];
    dtlims = [dtlims,dtlim'];
end

% Right now D1 is a stack of time series. The next step is to construct a
% whole-domain vector that varies first in space, then in time. A
% complication is the fact that the records each have different lengthts in
% time.
[Dmat,cot] = yv2cores(D1,dtlims,TSTEP);

[m n] = size(Dmat);

% create indices that relate the different records to the single,
% whole-domain vector. This looks tricky, but it works!
% NB: D1 - Dwd(yi(yi~0)) = 0
a = ones(size(Dmat));
a(isnan(Dmat)) = 0;
at = a';
yi = reshape(cumsum(at(:)),n,m)';
yi(~a) = 0;

% create the inverse index set to yi:
yiv = yi(:);
lyiv = length(yiv(yiv~=0));
yi_inv = zeros(lyiv,1);
cv = 1:lyiv;
for ii = 1:lyiv
    if any(yiv==ii)
        yi_inv(ii) = cv(yiv(yiv~=0)==ii);
    else
        yi_inv(ii) = 0;
    end
end

% now create the whole-domain vector Dwd and associated time vector Dwdt
Dmatt = Dmat';
bv = isnan(Dmatt(:)); % logical vector showing where data are missing. important for adjusting the design matrix!
Dwd = Dmatt(:);
Dwd(bv) = [];
Dwdtm = repmat(cot,n,1);
Dwdt = Dwdtm(:);
Dwdt(bv) = [];

% convert the matrix P to be ordered by space, then time
Pm = Pm(:,yi_inv);
Pm = Pm(yi_inv,:);

% convert yi to logical form yil
yil = zeros(length(Dwd),n);
for ii = 1:n
    yii = yi(:,ii);
    yii(~yii) = [];
    yil(yii,ii) = 1;
    yil = ~~yil;
end

LGMlims = [-25000 -20000];
tLGM = (tp>LGMlims(1) & tp<LGMlims(2));
LGMmv = nanmean(Dmat(tLGM,:)) + d18mv';

save Dlsfile Dwd Dwdt tp yil dtlims core_order bv Pm d18mv LGMmv 
