% tvsvd.m and plot_tvsvd.m make the incorrect assumption that the scaled resolution matrices are symmetric. Here we correct that assumption.
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
%[uv sv vv] = svds(reshape((bsxfun(@times,v230,cSd(:))),2806,[]),10);
%[uv sv vv] = svds(reshape(sparse(bsxfun(@times,v230,cSd(:))),2806,[]),10);

[q r] = qr(bsxfun(@times,v230,cSd(:)),0);
cSdi = 1./cSd(:); cSdi(cSdi==inf) = 0;
[p t] = qr(bsxfun(@times,v230,cSdi),0);
[uq sq vq] = svd(r*t');

% p*vq are the right singular vectors of Tv. Here we find their dominant spatial patterns
[urtv srtv vrtv] = svds(reshape(sparse((chol(sq)*(p*vq)')'),2806,[]),10);
% q*uq are the left singular vectors of Tv. Here we find their dominant spatial patterns
[ultv sltv vltv] = svds(reshape(sparse(q*uq),2806,[]),10);

%% Now project a pattern onto Tv
resp = q*uq*chol(sq)*reshape(reshape(sparse(p*vq),2806,[])'*urtv(:,1),230,101);
% is it dominated by a single spatial pattern? i.e. does sresp have one
% value much larger than all the others?
[uresp sresp vresp] = svds(reshape(resp,2806,[]),10);

%% Now project a pattern onto Tv transpose
resp2 = reshape(reshape(sparse(q*uq),2806,[])'*ultv(:,2),230,101)'*vq'*p';
% is it dominated by a single spatial pattern?
[uresp2 sresp2 vresp2] = svds(reshape(resp2,2806,[]),10);

%% Project an interhemispheric pattern
load ../../tmi/tracerobs_4deg_33lev_woce
ip = plotsurf(2*((LAT(jt)>0)' & kt==1)-1);
resp3 = q*uq*reshape(reshape(sparse(p*vq),2806,[])'*ip(~isnan(ip)),230,101);
% is it dominated by a single spatial pattern? i.e. does sresp have one
% value much larger than all the others?
[uresp3 sresp3 vresp3] = svds(reshape(resp3,2806,[]),10);

%% Trying something else
resp4 = q*uq*reshape(reshape(sparse(p*vq),2806,[])'*urtv(:,1),230,101);
[uresp4 sresp4 vresp4] = svds(reshape(resp4,2806,[]),10);

%% Freq dependence?

%% which patterns are best resolved?





