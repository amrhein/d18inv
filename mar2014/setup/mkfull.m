function gfulls = mkfull

TSTEP = 200;
ddir = 'cadjs/';
adjfiles = dir([ddir,'*.mat']);
lc = numel(adjfiles); % number of obs
load ../../tmi/tracerobs_4deg_33lev_woce

% create a matrix for mixing surface conditions through the ML
Mx = nan(length(it),2806);
ns = find(kt==1);
nsn = sum(kt==1);
for jj = 1:2806
    Mx(:,jj) = (it==it(ns(jj)) & jt==jt(ns(jj)) & inmixlyr);
end
MxMl = Mx(inmixlyr,:);

r1 = [];
for ii = 1:lc
    
    load([ddir,adjfiles(ii).name])
    
    Ti = T(1):TSTEP:T(end);
    
%    Cs = C(:,inmixlyr); % done already
    Csd = DGradient(C,T,1,'2ndOrder');
    Csd(1,:) = 0; % require that inl times have 0 conc
    Csd2 = DGradient(Csd,T,1,'2ndOrder');
    Csd2(1,:) = 0;% require that inl times have 0 conc
    Csd2i = interp1(T,Csd2,Ti);
    nc = cumsum(Csd2i);
    nm = sum(nc(end,:));
    
    % the normalization step
    Csd2in = Csd2i/nm;
    
    % sum over the vertical
    Csd2inv = (Csd2in*MxMl)';
    % pull out one surface state after another. ordered by space first, then time
    r1 = [r1,Csd2inv(:)];
    
end

m = length(r1)/nsn; % number of time steps

% make gfull

gfull = nan(m,lc*nsn);

for ii = 1:m
    inds = (1:nsn)+(ii-1)*nsn;
    gfull(ii,:) = reshape(r1(inds,:),[],1)';
end

% sparsify

gfulls = gfull;
gfulls(gfulls<10e-10)=0;
gfulls = sparse(gfulls);

save gfulls gfulls

