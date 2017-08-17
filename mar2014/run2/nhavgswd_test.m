% nhavgswd
% this is a whale of a script

%% Setup
clear
close all

k = 230; % k' used in ssvd
km = 99; % k' used in the overdetermined problem of estimating NH
          % and SH means (or really, over any region as specified below)

addpath ../../../export_fig

fl = ['solout',num2str(k)];
load(fl)

load ../../tmi/c_all_4deg.mat
load ../../tmi/tracerobs_4deg_33lev_woce.mat

% N = 25;
xestmat = reshape(xest,2806,[]);

[lonmg,latmg] = meshgrid(2:4:358,-90:4:88);

addpath ../utils/
addpath ../analysis/
load ../setup/Dlsfile

ld = length(Dwd);
nc = numel(core_order);

% Cholesky decompose the error covariance matrix
cPm = chol(Pm);

load ../setup/gfulls
stlims = [-25000 -5000];
TSTEP = 200;
nt = 101;
% obtain the whole-domain design matrix
nr = length(gfulls)/nc;
 Gf = mkG(stlims,TSTEP,gfulls,nc,nr);
 
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
%% Define some regions. In reginds, the northern box is in column 1
%% and the southern in column 2.

indmap = plotsurf(0*var(xestmat')); close
indmap(~isnan(indmap)) = 1:2806;
reginds = false(2806,1);

% Define the Northern box
nb = indmap(~isnan(indmap) & (latmg>=46 & latmg<50) & ((lonmg>=0 & lonmg<=30) |...
      (lonmg>=290 & lonmg<=360)));
 reginds(nb,1 ) = 1;
 plotsurf(reginds(:,1));

% Define the Southern box
sb = indmap(~isnan(indmap) & (latmg>=-84 & latmg<-60) & ((lonmg>=0 & lonmg<=30) |...
      (lonmg>=290 & lonmg<=360)));
 reginds(sb,2) = 1;
 plotsurf(reginds(:,2));

% index the regions in the new dimension
nsinds = reginds;
rib = logical(sum(reginds,2)); %places/times for either region
nsinds(~rib,:) = [];
nns = sum(reginds(:));
nn = sum(reginds(:,1));
ns = sum(reginds(:,2));

%% Derive the solution covariance submatrices for these two regions,
% including their cross-covariance, at every reconstruction time
%if ~exist('Mr','var')
mfn = ['Mout' num2str(k)];    
load(mfn)
%end
%%
[r c] = size(M);
Mr = reshape(bsxfun(@times,M,cSd'),2806,[]);
%Mr = reshape(M,2806,[]);
%clear M

% Derive the covariance matrix factors for the two regions...
Mrns = Mr(~~sum(reginds,2),:);
Mns = reshape(Mrns,[],c);
% find the svd of Mns*Mns' to generate its eigenvectors
[q,r] = qr(Mns,0);
[um sm ~] = svd(r*r');
e = q*um; % eigenvectors
er = e(:,1:k); % others are 0 within eps
                 %Wic = diag(diag(sm(1:k,1:k)).^-1/2)*uw(:,1:k)';

% Construct the design matrix
md = sparse(nns*nt,2*nt);
for ii = 1:nt
    indr = (ii-1)*nns+1;
    indc = (ii-1)*2+1;
    md(indr:(indr+nn-1),indc) = 1;
    md((indr+nn):(indr+nns-1),indc+1) = 1;
end

% These are the elements of xest that fall in both regions; they
% are the 'data' in this case. 
xestns = reshape(xestmat(rib,:),[],1);
  
% Project the eigenvectors of the covariance matrix onto the
% problem
y = er'*xestns;
E = er'*md;

% should I have column normalized here?
% I need to column normalize so that the region with more points
% doesn't 'win'...
ESd = 1./sqrt(sum(E.^2,1));
% diagonal of cholesky decomp of S
EcSd = sqrt(ESd);
EcSd(isnan(EcSd) | EcSd==inf) = 0;
% diagonal of inverse cholesky decomp of S
EcSdi = sqrt(sqrt(sum(E.^2,1)));
EcSdi(EcSdi==inf) = 0;

Es = E*diag(EcSd);

% now the problem is just y = Em + phi. <phi phi'> is
% diagonal and the diagonal elements are the eigenvalues of <theta
% theta'>

% row weighting matrix
Wic = diag(diag(sm(1:k,1:k)^(-1/2)));
wE = Wic*E;

%%
[uM sM vM] = svd(Wic*E);
uMk = uM(:,1:km);
sMk = sM(1:km,1:km);
vMk = vM(:,1:km);
mestu = tsvd(uM,sM,vM,inv(Wic)*y,km);
mest = diag(EcSd)*mestu; % undoing column scaling
nh = mest(1:2:end);
sh = mest(2:2:end);
time = -25000:200:-5000;

% soln covariance
Cmm = diag(EcSd)*vMk*inv(sMk)*uMk'*Wic*Wic*uMk*inv(sMk)*vMk'*diag(EcSd)';
% standard errors
dcmm = diag(Cmm);
nhse = sqrt(dcmm(1:2:end));
shse = sqrt(dcmm(2:2:end));

figure(2)
clf
hold all
errorbar(time,nh,nhse)
errorbar(time,sh,shse)
legend('nh','sh')

% Now for the difference of the problem.
Dt = sparse(nt*2,nt);
for ii = 1:nt
    rint = (ii-1)*2+1;
    Dt(rint,ii) = 1;
    Dt(rint+1,ii) = -1;
end
D = Dt'; % whoops

d = D*mest;
Cdd = D*Cmm*D';
% enough of this. let's just chop off the last two times for now
% (they are zero - in the nullspace of the problem upstream).
dm = d(1:end-2);
Cddm = Cdd(1:end-2,1:end-2);
Cddc = chol(Cddm);

% Now the problem is just-determined and the unweighted design
% matrix is just identity. After weighting, it is inv(Cyyc).
%
[ud sd vd] = svd(inv(Cddc));
udk = ud(:,1:kd);
sdk = sd(1:kd,1:kd);
vdk = vd(:,1:kd);
dest = tsvd(ud,sd,vd,inv(Cddc)*dm,kd);
%dest = tsvd(ud,sd,vd,dm,70);
%plot(dest)
Cdede = vdk*inv(sdk)*udk'*Cddm*udk*inv(sdk)*vdk';
dese = sqrt(diag(Cdede));

figure(3)
clf
errorbar([-25000:200:-5400],dest,dese)
ylabel('NH is colder/saltier -->')



