% nhavgs

%% Setup
clear
close all

k = 230;
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
nb = indmap(~isnan(indmap) & (latmg>=46 & latmg<80) & ((lonmg>=0 & lonmg<=30) |...
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

%% Derive the solution covariance submatrices for these two regions,
% including their cross-covariance, at every reconstruction time
if ~exist('Mr','var')
    load Mout230
end
%%
[r c] = size(M);
Mr = reshape(bsxfun(@times,M,cSd'),2806,[]);
%Mr = reshape(M,2806,[]);
%clear M

% Derive the covariance matrix factors for the two regions...
Mrns = Mr(~~sum(reginds,2),:);
% I'm interested in the difference!
Mrns(reginds(2,:),:) = -Mrns(reginds(2,:),:)/sum(reginds(:,2));
%Mrns(reginds(2,:),:) = -Mrns(reginds(2,:),:);
Mrns(reginds(1,:),:) = Mrns(reginds(1,:),:)/sum(reginds(:,1));
Mns = reshape(Mrns,[],c);
% find the eigenvalues of Mns*Mns'
[q,r] = qr(Mns,0);
[um sm ~] = svd(r*r');
e = q*um; % eigenvectors
e230 = e(:,1:230); % others are 0 within eps

% now the problem is just y = Dm + theta
xestns = reshape(xestmat(rib,:),[],1);
timp = reshape(mean(reshape(e230,sum(rib),[])),[],230);
y = e230'*xestns;

% this is just the unweighted average!!
plot(-25000:200:-5000,timp*y)


% ...for the NH...
Mrn = Mr(reginds(:,1),:);
Mn = reshape(Mrn,[],c);

% ...and for the SH.
Mrs = Mr(reginds(:,2),:);
Ms = reshape(Mrs,[],c);

%% Construct the averaging operators for the two regions and a difference
% operator between the two regions, at each time
nt = 101;
nns = sum(reginds(:));
nn = sum(reginds(:,1));
ns = sum(reginds(:,2));

% initialize operators. Each one is just a column; see 2.7.5 in
% Wunsch 2006
mns = [];
mn = [];
ms = [];
% initialize uncertainty (one number for each time)
uns = [];
un = [];
us = [];

%    ep = 1e-2; % pinv tolerance

    % Derive the differencing operator at each time for each of the
    % three desired quantities (NH, SH, and their difference)

    % NH operator
    %    Mnt = Mn((ii-1)*nn+1:ii*nn,:);
    %D = ones(nn,1);
    %mn(:,ii) = pinv(D'*pinv(Mnt*Mnt',ep)*D,ep)*D'*pinv(Mnt*Mnt',ep);
    %un(ii) = 1./sum(sum(pinv(Mnt*Mnt',ep)));
    
    % SH operator
    % Mst = Ms((ii-1)*ns+1:ii*ns,:);
    % D = ones(ns,1);
    % ms(:,ii) = pinv(D'*pinv(Mst*Mst',ep)*D,ep)*D'*pinv(Mst*Mst',ep);
    % us(ii) = 1./sum(sum(pinv(Mst*Mst',ep)));

    % differencing operator
    %     D = nsinds(:,1)-nsinds(:,2);

    %     Mnst = Mns((ii-1)*nns+1:ii*nns,:);
    % mns(:,ii) = pinv(D'*pinv(Mnst*Mnst',ep)*D,ep)*D'*pinv(Mnst*Mnst',ep);
    % uns(ii) = 1/(D'*pinv(Mnst*Mnst',ep)*D); % not the sum!


% Apply the operator to the solution to obtain averages and uncertainty
% estimates

% north
na = sum(xestmat(reginds(:,1),:).*mn);

% south
sa = sum(xestmat(reginds(:,2),:).*ms);

% north-south
nsd = sum(xestmat(~~sum(reginds,2),:).*mns);

close
t = -25000:200:-5000;
%plot(t,na,t,sa,t,nsd)

close
hold all
errorbar(t,na,sqrt(un))
errorbar(t,sa,sqrt(us))
errorbar(t,nsd,sqrt(uns))
axis tight
legend('SH','NH','NH-SH')
title('Covariance weighted average')
export_fig('-pdf',['Figs/nhavg',num2str(log10(ep))])

%subplot(2,1,2),plot(t,nanw,t,sanw,t,nsdnw)
%legend('SH','NH','NH-SH')
%title('Unweighted average')

