% Solves the whole #!
% ssvd_wm.m: also solves for the mean of the problem.

%% Setup

addpath ../utils/

clear

k = [300] % svd truncation parameter

TSTEP = 200;
stlims = [-25000 -5000];
nc = 8; % number of sediment cores

% time axis for the reconstruction
t = stlims(1):TSTEP:stlims(2);
nt = length(t);

ntg = 25; % number of time steps in the gn fn vectors

fsp = [1075 1439 1550 821]; % full screen position; for plotting

% Obtain sediment core records in whole-domain form
load ../setup/Dlsfile

ld = length(Dwd);
M = numel(core_order);

load ../setup/gfulls
nr = length(gfulls)/M; % number of surface grid boxes

% Cholesky decompose the error covariance matrix
cPm = chol(Pm);

%% Cook up Gbar

% obtain the whole-domain design matrix
Gf = mkG(stlims,TSTEP,gfulls,M,nr);
 
% % correct for different record lengths
Gb = Gf(~bv,:);

% weight by Pm
Gbp = sparse(cPm\Gb);
 
% use this to unweight the (column scaling in the) solution
% diagonal of column scaling matrix S
Sd = 1./sqrt(sum(Gbp.^2,1));
% diagonal of cholesky decomp of S
cSd = sqrt(Sd);
cSd(isnan(cSd) | cSd==inf) = 0;

%load Gbp
load('Gbpusv_sqrt.mat');

% subselect for the IC problem
G0 = Gf(:,1:nr*50);

% create the matrix for the IC problem
% Gf does not have any of the scaling applied in Gbpusv_sqrt
Gicb = Gf(end-8:end,:);%1:nr*50);
Gic = nan(M,nr);
for ii = 1:M
    Gic(ii,:) = sum(reshape(Gicb(ii,:),nr,[]),2);
end

% row scaling matrix W
sig = 0.2; 
Wd0 = diag(ones(1,nc) * sig); % so long as this is diagonal, row weighting does not matter
% its Chol
Wdoc = diag(diag(Wd0).^(1/2));
% its inverse Cholesky decomposition
Wdoci = diag(diag(Wd0).^(-1/2));
Gicw = Wdoci*Gic;

% diagonal of column scaling matrix S
iSd = 1./sqrt(sum(Gicw.^2,1));
% diagonal of cholesky decomp of S
ciSd = sqrt(iSd);
Gicwc = bsxfun(@times,Gicw,ciSd);
Gicwc(isnan(Gicwc)) = 0;

%% solve for ICs
kic = 4;
% put the mean in here!
ymat = yv2coreswd(Dwd,yil);
%yold = mean(ymat(1:10,:))';
yold = mean(ymat(1:10,:))';
yhol = nanmean(ymat((end-10):end,:))';
[uic,sic,vic] = svd(Gicwc);
yolds = Wdoci*yold;
[xicu] = tsvd(uic,sic,vic,yolds,kic);
% unscaling
xic = ciSd'.*xicu;
% remember: cSd has some infinite values!
xicn = xic; xicn(isnan(xicn)) = 0;
yic = Wdoc*Gicwc*xicu;

% find the soln covariance
%Mic = vic(:,1:kic)*(sic(1:kic,1:kic)\uic(:,1:kic)')*Wdoc';
Mic = vic(:,1:kic)*(sic(1:kic,1:kic)\uic(:,1:kic)');
Cxxicu = Mic*Mic';
% undo column scaling
Cxxic = diag(ciSd)*Cxxicu*diag(ciSd);
Cxxic(isnan(Cxxic)) = 0;

% generate vectors of forward run output at the core locations under a
% forcing of 5kyr initial conditions followed by 5kyr of nothing
xic10k = G0*[repmat(xic,25,1);zeros(nr*25,1)];

hl = 25*M;

xic10k2h = [xic10k(hl+1:end);zeros(hl,1)];

fname = ['icout' num2str(kic)]
%save(fname, 'xic10k2h', 'xic', 'yic', 'yold','sic','uic','vic','kic','Cxxic','ciSd')

%% solve for the mean problem

% create the matrix for the mean problem
Gm = reshape(sum(gfulls),nr,nc)';

% don't weight the problem by the mean uncertainties
Gmp = Gm;%chol(diag(d18mP))\Gm;

% weight the columns of the problem
% diagonal of column scaling matrix S
mSd = 1./sqrt(sum(Gmp.^2,1));
% diagonal of cholesky decomp of S
cmSd = sqrt(mSd);

% apply column weighting
Gmpcw = bsxfun(@times,Gmp,cmSd);
Gmpcw(isnan(Gmpcw)) = 0;

km = 4;
[um,sm,vm] = svd(full(Gmpcw));
[xmu] = tsvd(um,sm,vm,d18mv,km);

% unscaling
xmn = cmSd'.*xmu;
% remember: cSd has some infinite values!
xm = xmn; xm(isnan(xm)) = 0;
yme = Gmpcw*xmu;

%save xm yme d18m km sm


%% solve the time-varying problem
k = 230;
% no mean:
%yp = cPm\(Dwd-xic10k2h(~bv));

% put the mean back in:
yp =  cPm\(Dwd-xic10k2h(~bv)+sum(bsxfun(@times,d18mv',yil),2));

% put the mean back in but subtract c:
c = 3.5
yp =  cPm\(Dwd-xic10k2h(~bv)+sum(bsxfun(@times,d18mv',yil),2)-c);

% only the mean:
%yp =  cPm\(sum(bsxfun(@times,d18mv',yil),2));

% constant values everywhere:
%yp =  cPm\(sum(bsxfun(@times,d18mv',yil),2));
%yp = cPm\(yp*0+1);

[xest] = tsvd(u,s,v,yp,k);


% unscaling
xest = cSd'.*xest;
% remember: cSd has some infinite values!
xest(isnan(xest)) = 0;

yest = Gbp*xest;
yestu = yest;
yestu = cPm*yest + xic10k2h(~bv);
%yestu = cPm*yest;

fname = ['solout_wm',num2str(k)]
save(fname,'yp', 'xest', 'yest', 'yestu', 's')

%
figure(1)
clf
% with the mean
%for ii= 1:8, subplot(3,3,ii),hold all,  plot(yestu(yil(:,ii))),
%plot(Dwd(yil(:,ii))+d18mv(ii)), end

% with the mean, subtracting a constant
for ii= 1:8, subplot(3,3,ii),hold all,  plot(yestu(yil(:,ii))+c),
 plot(Dwd(yil(:,ii))+d18mv(ii)), end

% without the mean
%for ii= 1:8, subplot(3,3,ii),hold all,  plot(yestu(yil(:,ii))), plot(Dwd(yil(:,ii))), end

% only the mean
%moy = sum(bsxfun(@times,d18mv',yil),2);
%for ii= 1:8, subplot(3,3,ii),hold all,  plot(yestu(yil(:,ii))), 
%plot(moy(yil(:,ii)))
%end

% constant values
%for ii= 1:8, subplot(3,3,ii),hold all,  plot(yestu(yil(:,ii))), 
%end

addpath ../analysis
%plot_xest_ts_novb_wm(k)

xestm = reshape(xest,2806,[]);
%figure(2),close,figure(2)
%for ii = 1:101, plotsurf(xestm(:,ii));caxis([-15 15]), colorbar, pause(0.2), end
