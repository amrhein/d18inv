% plot_tvsvd.m


%%
clear
%% Setup
%matlabpool

addpath ../utils/
addpath ../analysis/
load ../setup/Dlsfile

ld = length(Dwd);
M = numel(core_order);

% Cholesky decompose the error covariance matrix
cPm = chol(Pm);

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
cSdi = sqrt(sqrt(sum(Gbp.^2,1)));
cSdi(cSdi==inf) = 0;

%load Gbp
load('Gbpusv_sqrt.mat');

v230 = v(:,1:230);

%% plot
close all
set(gcf,'color','w','position',[1 1 1024 1201])

quqr = reshape(bsxfun(@times,v230,cSd(:)),2806,[]);
pvqr = reshape(bsxfun(@times,v230,cSdi(:)),2806,[]);
%% Now: pointwise
gainmat = [];
% parfor ii = 1:2806
for ii = 1:2806    
    Us = reshape(quqr(ii,:),101,230);
    Vs = reshape(pvqr(ii,:),101,230);
    tf = Us*Vs';
    keyboard  
    % plot the temporal filter gain
    [f,gain] = getgain(tf(50,:),1/200);
    gainmat(:,ii) = gain;
    plot(f,gain)
end

%save gainmat f gainmat

bubbsurf(gainmat(1,:),10e-4,300);
export_fig('pdf','Figs/gain0')
