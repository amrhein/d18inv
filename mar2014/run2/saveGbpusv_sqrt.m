% Modified from saveGbpusv DEA 26 Mar 2014 to include column scaling as described in Wunsch 2006 p 108-110

%% Setup
clear

addpath ../setup

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
load Dlsfile

ld = length(Dwd);
M = numel(core_order);

% Cholesky decompose the error covariance matrix
cPm = chol(Pm);

%% Cook up Gbar

%load ../23july/4weddell/gfulls23July
load gfulls

% obtain the whole-domain design matrix
nr = length(gfulls)/M;
 Gf = mkG(stlims,TSTEP,gfulls,M,nr);

% % correct for different record lengths
 Gb = Gf(~bv,:);

% % weight by Pm
 Gbp = sparse(cPm\Gb);

% weight by sqrt of column length
%S = sparse(diag(1./sqrt(sum(Gbp.^2,1))));
%Gbpc = Gbpc*chol(S);
%lrt = full(sum(Gbp.^2).^(1/4));
% diagonal of column scaling matrix S
Sd = 1./sqrt(sum(Gbp.^2,1));
% diagonal of cholesky decomp of S
cSd = sqrt(Sd);
Gbpc = bsxfun(@times,Gbp,cSd);
Gbpc(isnan(Gbpc)) = 0;

%% expensive! only do this once!
   [u s v] = svd(full(Gbpc),'econ');
   %vs = sparse(v);
%save('-v7.3','Gbpusv.mat','u','s','vs');
save('-v7.3','Gbpusv_sqrt.mat','u','s','v');

