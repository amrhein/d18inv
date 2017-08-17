function [co,cot] = yv2cores(yv,dtlims,TSTEP)
%function [co,cot] = yv2cores(yv,dtlims,TSTEP)
% A function that takes the time ranges of core records, the timestep,
% and the whole domain vector y and returns a matrix of nans and
% core records all on the same time axis as well as the accompanying time
% vector

% compute the starting and ending indices for each d18O record
lens = diff(dtlims)/TSTEP;
ts = lens+1;
tsc = cumsum(ts);
starts = [1,tsc(1:end-1)+1];
ends = tsc;

% whole time duration
cot = min(dtlims(1,:)):TSTEP:max(dtlims(2,:));
nt = length(cot);
co = nan(nt,length(dtlims)); %initialize matrix of cores

% populate the matrix co
for i = 1:length(dtlims)
    inds = starts(i):ends(i);
    tdt = (dtlims(1,i):TSTEP:dtlims(2,i));
    co(cot>=min(tdt) & cot<=max(tdt),i) = yv(inds);
end