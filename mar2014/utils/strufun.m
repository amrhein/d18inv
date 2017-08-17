function [blag,brms,bvar,lagv,rmsv,lagm] = strufun(t,s,N,nb)
% Computes the structure function for a time series.
% From Press et al.1992, 'Time Delay of Gravitational Lens 0957 + 561 I.
% Eq. 22
% t is a vector of times
% s is a vector of measured values
% e is a vector of measured value uncertainties
% nb is the number of bins used

l = length(t);
rmsm = [];
lagm = [];

for ii = 2:l
    for jj = 1:ii
        lagm(ii,jj) = abs(t(ii)-t(jj));
        rmsm(ii,jj) = (s(ii)-s(jj))^2 - N(ii,ii) - N(jj,jj);
    end
end

rmsm(rmsm<0) = nan;
li = ~triu(ones(size(rmsm)));

rmsv = rmsm(li(:));
lagv = lagm(li(:));

b_e = logspace(log10(min(lagv)),log10(max(lagv)),nb+1); % bin edges
if any(lagv==0)
    error('no repeat measurements please!')
end
brms = nan(nb,1); % binned avg rms
bvar = brms; % intra-bin variance of binned rms
blag = b_e(1:end-1)+diff(b_e)/2; % bin centers

for ii = 1:nb
    brms(ii) = nanmean(rmsv(lagv>=b_e(ii) & lagv<b_e(ii+1)));
    bvar(ii) = nanvar(rmsv(lagv>b_e(ii) & lagv<b_e(ii+1)));
end
    
brms = brms(:);
blag = blag(:);
bvar = bvar(:);

if all(isnan(blag))
    keyboard
end
    
    
    
    

