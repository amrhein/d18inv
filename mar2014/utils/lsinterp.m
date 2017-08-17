% lsinterp.m

function [yi,P,S,ym,yms] = lsinterp(y,t,tp,N)
% yi are values at the interpolant times tp
% P is the minimum uncertainty of the interpolated record
% S is the estimate of the signal (minus noise!) covariance computed by
% fiting an exponent to the structure function

[ym,yms] = get_ym(y,t,N);

ly = length(y);

ymm = y - ym;
msy = mean(ymm.^2);

% compute the structure function estimate
[blag,brms,bvar,lagv,rmsv,lagm] = strufun(t,ymm,N,20);
bg2 = (blag>0 & ~isnan(brms));
p = polyfit(log10(blag(bg2)),log10(brms(bg2)),1);
a = 10^p(2);
b = p(1);

%disp(['a = ' num2str(a) ', b = ' num2str(b)])

strf = @(tau) a*tau.^(b);

% lagm is a lower triangular matrix. make it full with a diagonal of zeros:
lagmf = lagm+lagm';
lagmf(~~eye(size(lagmf))) = 0;

% generate the signal covariance function
S = msy - 0.5*strf(abs(lagmf));
%S(S(:)<0)=0;

SN = (S+N);

% Generate the correlation matrix Ss (S* in Rybicki and Press 1992 eqn 5)
Ss = [];

for ii = 1:ly
    Ss(ii,:) = msy - 0.5*strf(abs(t(ii) - tp));
end
%Ss(Ss(:)<0)=0;


yi = Ss'*(SN\ymm);
yi = yi(:);

% Generate an estimate of the solution covariance matrix at the
% interpolated points
lyp = length(tp);
Rxx = [];
for ii = 1:lyp
    Rxx(ii,:) = msy - 0.5*strf(abs(tp(ii) - tp));
end

P = Rxx - Ss'*(SN\Ss); % minimum uncertainty




