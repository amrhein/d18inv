%% Setup

clear

load Vtot_DEA_30Sept.mat 

TSTEP = 200;
% if you change stlims, you will have to recompute Gbusv!
stlims = [-25000 -5000];
sig = 0.2; % permil std dev. this is also built into mkD.m

% time axis for the reconstruction
t = stlims(1):TSTEP:stlims(2);
nt = length(t);

ntg = 25; % number of time steps in the gn fn vectors
gam = 1; % tapering parameter

fsp = [1075 1439 1550 821]; % full screen position; for plotting

%% Obtain sediment core records in whole-domain form
load('../23july/4weddell/Dlsfile');

ld = length(Dwd);
M = numel(core_order);

% Cholesky decompose the error covariance matrix
cPm = chol(Pm);

%% Cook up Gbar

load ../23july/4weddell/gfulls23July

% obtain the whole-domain design matrix
nr = length(gfulls23July)/M;
 Gf = mkG(stlims,TSTEP,gfulls23July,M,nr);

% % correct for different record lengths
 Gb = Gf(~bv,:);

% % weight by Pm
 Gbp = sparse(cPm\Gb);

% keep the 600 leading surface points. They fill 95% of the ocean!
%Vtotv = Vtot(~isnan(Vtot))/nansum(Vtot(:));
Vtotv = Vtot(~isnan(Vtot));
[vsorted,sorti] = sort(Vtotv,'descend');
keepv = sorti(1:600);
tp = zeros(2806,1); tp(keepv)=1;
tpr = repmat(~~tp',1,101);
Gbpt = Gbp(:,tpr);
Gbptv = bsxfun(@times,Gbpt,repmat(1./Vtotv(~~tp)',1,101));

%Gbpv = bsxfun(@times,Gbp,repmat(1./Vtotv',1,101));


%% expensive! only do this once!
   [u s v] = svd(full(Gbptv),'econ');
%   vs = sparse(v);
save('-v7.3','4weddell/Gbpusv30Sept_V2.mat','u','s','v');

