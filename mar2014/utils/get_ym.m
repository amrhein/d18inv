% lsinterp.m

function [ym,yms] = get_ym(y,t,N)

msy = mean(y.^2);
ly = length(y);
ynm = y-nanmean(y);

% compute the structure function estimate
[blag,brms,bvar,lagv,rmsv,lagm] = strufun(t,ynm,N,20);
bg2 = (blag>=0 & ~isnan(brms));
p = polyfit(log10(blag(bg2)),log10(brms(bg2)),1);
a = 10^p(2);
b = p(1);

strf = @(tau) a*tau.^(b);

% lagm is a lower triangular matrix. make it full with a diagonal of zeros:
lagmf = lagm+lagm';
lagmf(~~eye(size(lagmf))) = 0;

% generate the signal covariance function
S = mean(ynm.^2) - 0.5*strf(abs(lagmf));
S(abs(S)==inf)=0; % can happen for blue spectra
%S(S(:)<0)=0;

% compute an oft-used inverse
iSN = inv(S+N);

% mean estimator. See Rybicki and Press 1992 or Wunsch 2006 2.413
ee = ones(ly,1);
ym = ee'*(iSN*y)/(ee'*iSN*ee);
yms = 1/(ee'*iSN*ee);


