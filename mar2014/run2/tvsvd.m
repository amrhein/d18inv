
%% Setup

addpath ../utils/

clear

load ../setup/Dlsfile

ld = length(Dwd);
M = numel(core_order);

% Cholesky decompose the error covariance matrix
cPm = chol(Pm);

%% Cook up Gbar

load ../setup/gfulls
stlims = [-25000 -5000];
TSTEP = 200;
% obtain the whole-domain design matrix
nr = length(gfulls)/M;
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

v230 = v(:,1:230);

% [uv sv vv] = svds(reshape(sparse(bsxfun(@times,v230,cSd(:))),2806,[]),10);

