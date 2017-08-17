% nhavgswd
% this is a whale of a script

%% Setup
clear
close all

k = 790; % k' used in ssvd
%km = 99; % k' used in the overdetermined problem of estimating NH
          % and SH means (or really, over any region as specified below)

nr = 2806; %number of surface regions
nto = 101; % number of reconstruction times
nt = 99; % number of times considered here. the last 2 are in the nullspace.
addpath ../../../export_fig

fl = ['solout',num2str(k)];
load(fl)

load ../../tmi/c_all_4deg.mat
load ../../tmi/tracerobs_4deg_33lev_woce.mat

xestmat = reshape(xest,2806,[]);
xestmat(:,nt+1:nto) = [];

[lonmg,latmg] = meshgrid(2:4:358,-90:4:88);

addpath ../utils/
addpath ../analysis/

%% Define some regions. In reginds, the northern box is in column 1
%% and the southern in column 2.

indmap = plotsurf(0*var(xestmat')); close
indmap(~isnan(indmap)) = 1:2806;
reginds = false(2806,1);

% Define the Northern box

% a small high-latitude box:
%nb = indmap(~isnan(indmap) & (latmg>=46 & latmg<80) & ((lonmg>=0 & lonmg<=30) | (lonmg>=290 & lonmg<=360)));

% the whole NH
 nb = indmap(~isnan(indmap) & latmg>0);

% a single location
%nb = indmap(~isnan(indmap) & latmg == 70 & lonmg==2);

reginds(nb,1 ) = 1;
plotsurf(reginds(:,1));

% Define the Southern box

% a small high-latitude box:
%sb = indmap(~isnan(indmap) & (latmg>=-84 & latmg<-60) & ((lonmg>=0 & lonmg<=30) | (lonmg>=290 & lonmg<=360)));

% the whole SH
sb = indmap(~isnan(indmap) & latmg<0);

% a single location
%sb = indmap(~isnan(indmap) & latmg == -70 & lonmg==2);

reginds(sb,2) = 1;
plotsurf(reginds(:,2));

% index the regions in the new dimension
nsinds = reginds;
rib = logical(sum(reginds,2)); %places/times for either region
nsinds(~rib,:) = [];
nns = sum(reginds(:));
nn = sum(reginds(:,1));
ns = sum(reginds(:,2));

% Define the area-weighting matrix

lat2806w =  abs(cos(pi/180*latmg(~isnan(indmap))));
% the diagonal:
%awd =[lat2806w(reginds(:,1));lat2806w(reginds(:,2))];
awd =[lat2806w(nb));lat2806w(sb)];

%% Load M
load Gbpusv_sqrt
Ma = v(:,1:k)*(s(1:k,1:k)\u(:,1:k)');
% eliminate the last two points
M = Ma;
M((end-2*nr+1):end,:) = [];
clear Ma;

[r c] = size(M);

% load the problem column scaling matrix
load Sd
% eliminate the last two points
Sda = Sd;
Sda((end-2*nr+1):end) = [];
Sd = Sda; clear Sda;
cSda = cSd;
cSda((end-2*nr+1):end) = [];
cSd = cSda; clear cSda;
ciSd = 1./cSd;
ciSd(ciSd==inf) = 0;

% Rescale M from the solution of the last problem
Mr = reshape(bsxfun(@times,M,cSd'),2806,[]);

%%

% Derive the covariance matrix factors for the two regions...
Mrns = Mr(rib,:);
Mns = reshape(Mrns,[],c);
% find the svd of Mns*Mns' to generate its eigenvectors
[q,r] = qr(Mns,0);
[um sm ~] = svd(r*r');
e = q*um; % left singular vectors
ke = min(k,min(size(e)));
er = e(:,1:ke); % others are 0 within eps

% Construct the design matrix
md = sparse(nns*nt,2*nt);
for ii = 1:nt
    indr = (ii-1)*nns+1;
    indc = (ii-1)*2+1;
    md(indr:(indr+nn-1),indc) = 1;
    md((indr+nn):(indr+nns-1),indc+1) = 1;
end

% I need to column normalize so that the region with more points
% doesn't 'win'...
mSd = 1./sqrt(sum(md.^2,1));
% diagonal of cholesky decomp of S
mcSd = sqrt(mSd);
mcSd(isnan(mcSd) | mcSd==inf) = 0;
% diagonal of inverse cholesky decomp of S
%mcSdi = sqrt(sqrt(sum(md.^2,1)));
%mcSdi(mcSdi==inf) = 0;

% Finally, scale the design matrix
mds = sparse(md*diag(mcSd));

% These are the elements of xest that fall in both regions; they
% are the 'data' in this case. 
xestns = reshape(xestmat(rib,:),[],1);
  
% Project the eigenvectors of the covariance matrix onto the
% problem
y = er'*xestns;
E = er'*mds;

% row weighting matrix given by the inverse cholesky decomposition
% of the singular values sm
W = sm(1:ke,1:ke);
Wic = chol(inv(W));

%
km = 198;
[uM sM vM] = svd(Wic*E);
uMk = uM(:,1:km);
sMk = sM(1:km,1:km);
vMk = vM(:,1:km);
mestu = tsvd(uM,sM,vM,Wic*y,km);
mest = diag(mcSd)*mestu; % undoing column scaling
nh = mest(1:2:end);
sh = mest(2:2:end);
time = -25000:200:-5000;

% soln covariance
%Cmm = diag(mcSd)*vMk*inv(sMk)*uMk'*W*uMk*inv(sMk)*vMk'*diag(mcSd)';
Cmm = diag(mcSd)*vMk*inv(sMk)*uMk'*uMk*inv(sMk)*vMk'*diag(mcSd)';
% standard errors
dcmm = diag(Cmm);
nhse = sqrt(dcmm(1:2:end));
shse = sqrt(dcmm(2:2:end));

figure(2)
clf
hold all
tis = [-25000:200:-5400]; % shorter time
errorbar(tis,nh,nhse)
errorbar(tis,sh,shse)
legend('nh','sh')

% Now for the difference of the problem.
Dt = sparse(nt*2,nt);
for ii = 1:nt
    rint = (ii-1)*2+1;
    Dt(rint,ii) = 1;
    Dt(rint+1,ii) = -1;
end
D = Dt';

d = D*mest;
Cdd = D*Cmm*D';
Cddc = chol(Cdd);

% Now the problem is just-determined and the unweighted design
% matrix is just identity. After weighting, it is inv(Cyyc).
%
kd = 99;
[ud sd vd] = svd(inv(Cddc));
udk = ud(:,1:kd);
sdk = sd(1:kd,1:kd);
vdk = vd(:,1:kd);
dest = tsvd(ud,sd,vd,inv(Cddc)*d,kd);
%Cdede = vdk*inv(sdk)*udk'*Cdd*udk*inv(sdk)*vdk';
% no column scaling to undo here
Cdede = vdk*inv(sdk)*udk'*udk*inv(sdk)*vdk';
dese = sqrt(diag(Cdede));
export_fig('-pdf',['Figs/nhavgs_whole_hemis_',num2str(k)])
set(gcf,'color','w')
%export_fig('-pdf',['Figs/nhavgs_highlat_',num2str(k)])

figure(1)
clf
errorbar([-25000:200:-5400],dest,dese)
ylabel('NH is colder/saltier -->')
set(gcf,'color','w')
%export_fig('-pdf',['Figs/nhdiffs_highlat_',num2str(k)])
export_fig('-pdf',['Figs/nhavgs_whole_hemis_',num2str(k)])


%% Exercise as a sanity check
% When the two regions are single points, the inferred nh and sh
% means should just give the same values and uncertainties as for
% the original two points! assuming the regions have been thus
% assigned, compare the output in figure 2 to the following:

% nh1 = xestns(1:2:end);
%sh1 = xestns(2:2:end);
%Cns = Mns*Mns';
%Cnsd = diag(Cns);
%nhse1 = sqrt(Cnsd(1:2:end));
%shse1 = sqrt(Cnsd(2:2:end));
%clf
%hold on
%errorbar(nh1,nhse1)
%errorbar(sh1,shse1)

%% Resolution!
lm = length(mds);
%xest = cSd'.*xest;
%mest = diag(mcSd)*mestu; % undoing column scaling
% the mean operator:
mo = diag(mcSd)*vMk*(sMk\uMk')*Wic*er';
mdsv = mo*(bsxfun(@times,cSd',v(1:lm,1:k)));
mdsvv = mdsv*(bsxfun(@times,ciSd',v(1:lm,1:k)))';
%mdsv = mo*(bsxfun(@times,1,v(1:lm,1:k)));
%mdsvv = mdsv*(bsxfun(@times,1,v(1:lm,1:k)))';
mdsvvmm = (md*diag(1./sum(md)))'-mdsvv;
m2 = mdsvvmm*mdsvvmm';

% now for the case that is completely spatially covarying within
% the two regions. This can be efficiently computed and stored by
% factoring into two parts. kron tiles this matrix, and the
% covariance matrix (which we never need to compute) is mbf*mbf'.
mbf = sparse(logical(kron(eye(nt),[[ones(nn,1);zeros(ns,1)],[zeros(nn,1);ones(ns,1)]])));
mmbf = mdsvvmm*mbf;
mb2 = mmbf*mmbf';

% Now compute the time covariance matrix of the interhemispheric
% differences for the two prior covariances.

mb2D = D*mb2*D';
m2d = D*m2*D';

