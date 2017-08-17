% Solves the whole shebang

%% Setup

addpath ../utils/

clear

TSTEP = 200;
stlims = [-25000 -5000];
sig = 0.2; % permil std dev. this is also built into mkD.m

% time axis for the reconstruction
t = stlims(1):TSTEP:stlims(2);
nt = length(t);

ntg = 25; % number of time steps in the gn fn vectors
gam = 1; % tapering parameter

fsp = [1075 1439 1550 821]; % full screen position; for plotting

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

%load Gbp
load('../setup/Gbpusv.mat');

% subselect for the IC problem
G0 = Gf(:,1:nr*50);

% create the matrix for the IC problem
Gicb = Gf(end-8:end,:);%1:nr*50);
Gic = nan(M,nr);
for ii = 1:M
    Gic(ii,:) = sum(reshape(Gicb(ii,:),nr,[]),2);
end

%% solve for ICs
kic = 3;
ymat = yv2coreswd(Dwd,yil);
yold = mean(ymat(1:10,:))';
[uic,sic,vic] = svd(Gic);
[xic] = tsvd(uic,sic,vic,yold,kic);
yic = Gic*xic;

% generate vectors of forward run output at the core locations under a
% forcing of 5kyr initial conditions followed by 5kyr of nothing
xic10k = G0*[repmat(xic,25,1);zeros(nr*25,1)];

hl = 25*M;
xic10k2h = [xic10k(hl+1:end);zeros(hl,1)];

fname = ['icout' num2str(kic)]
save(fname, 'xic10k2h', 'xic', 'yic', 'yold','sic','uic','vic','kic')


%% solve
k = [413]
yp = cPm\(Dwd-xic10k2h(~bv));
[xest] = tsvd(u,s,v,yp,k);
yest = Gbp*xest;
yestu = cPm*yest + xic10k2h(~bv);

fname = ['solout',num2str(k)]
save(fname,'yp', 'xest', 'yest', 'yestu', 's')

% Solution covariance matrix
% From Wunsch (2006) Eq. 2.293:
% Cxx = vk*inv(s)*uk'*<n*n'>*uk*inv(s)*vk'
% Here, I will only compute the diagonal.

% recall that chol(A)'*chol(A) = A.

% Construct the LHS of eq. 2.293
M = v(:,1:k)*(s(1:k,1:k)\u(:,1:k)')*cPm';

% Now compute the diagonal of the covariance matrix
Cxx = sum(M.^2,2);

fname = ['Cxxout',num2str(k)]
save(fname,'Cxx')

%% Resolution matrices

% Tv is too expensive to calculate entirely, and there's 
% really no need to anyway. Instead, calculate a subset:
 
% inds = 1:nr;
% % index in time corresponding to the delta function 
% tvi = 50; 
% iii = tvi*nr+inds;
% 
% % number of time steps on either side of the impulse
% nes = 2;
% 
% % indices for the rows of Tv
% tvii = (2806*(tvi-nes)+1):(2806*(tvi+nes));
% 
% Tv = (v(tvii,1:k)*v(iii,1:k)');
% 
% Tu = (u(:,1:k)*u(:,1:k)');
% 
% fname = ['resout',num2str(k)]
% 
% save(fname,'-v7.3','Tv','Tu','k','nes','tvi')
