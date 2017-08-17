function [co] = yv2coreswd(yv,yil)
%function [co,cot] = yv2cores(yv,dtlims,TSTEP)
% A function that takes the time ranges of core records, the timestep,
% and the whole domain vector y and returns a matrix of nans and
% core records all on the same time axis as well as the accompanying time
% vector

[m n] = size(yil);
sy = sum(yil);
co = nan(max(sy),n);
%sofar = [0,cumsum(sum(yil))];

for ii = 1:n
    %co(1+sofar(ii):sofar(ii+1),ii) = yv(yil(:,ii));
    co(1:sy(ii),ii) = yv(yil(:,ii));
end
