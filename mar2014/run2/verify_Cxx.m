% checks the size of the standard error Cxx for the paper

%% Setup

addpath ../utils/

clear

k = [230] % svd truncation parameter
%GAM = 0.5; % tapering parameter

TSTEP = 200;
stlims = [-25000 -5000];

% time axis for the reconstruction
t = stlims(1):TSTEP:stlims(2);
nt = length(t);

ntg = 25; % number of time steps in the gn fn vectors

%% Obtain sediment core records in whole-domain form
load ../setup/Dlsfile

ld = length(Dwd);
M = numel(core_order);

% Cholesky decompose the error covariance matrix
cPm = chol(Pm);

%% Cook up Gbar

load ../setup/gfulls

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

% Solution covariance matrix
% From Wunsch (2006) Eq. 2.293:
% Cxx = vk*inv(s)*uk'*<n*n'>*uk*inv(s)*vk'
% Here, I will only compute the diagonal.

% recall that chol(A)'*chol(A) = A.

% Construct the LHS of eq. 2.293
%M = v(:,1:k)*(s(1:k,1:k)\u(:,1:k)')*cPm';
M = v(:,1:k)*(s(1:k,1:k)\u(:,1:k)');
% unweight by Sd
Mu = bsxfun(@times,M,cSd');

Cxxd = sum(Mu.^2,2);
