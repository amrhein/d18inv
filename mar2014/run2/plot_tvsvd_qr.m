% plot_tvsvd.m


%%
clear
L = 5; % number of patterns to plot
%% Setup

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

[q r] = qr(bsxfun(@times,v230,cSd(:)),0);
[p t] = qr(bsxfun(@times,v230,cSdi(:)),0);
[uq sq vq] = svd(r*t');

% p*vq are the right singular vectors of Tv. Here we find their dominant spatial patterns
[urtv srtv vrtv] = svds(reshape(sparse(((p*vq)')'),2806,[]),10);

%% plot
close all
set(gcf,'color','w','position',[1 1 1024 1201])

% isolate the temporal filter part
% Vk = tf*tp (temporal filter times toeplitz)
%  tf = v230\toeplitz([u(:,1)',zeros(2806*100,1)],[u(1,1);zeros(100,1)])

% Assume that the temporal filter is nearly stationary

for ii = 1:L
    pii = (ii-1)*3 + 1; % subplot index
    subplot(L,3,pii:(pii+1))
    plotsurf(urtv(:,ii));
    shading flat
    axis equal
    axis tight
    colorbar
    % compute the temporal filter
    ps = urtv(:,ii); % spatial pattern
    
    Us = reshape(ps(:)'*reshape(q*uq,2806,[]),101,230);
    Vs = reshape(ps(:)'*reshape(p*vq,2806,[]),101,230);
    tf = Us*sq*Vs';
  
    keyboard
    % plot the temporal filter gain
    subplot(L,3,pii+2)
    [f,gain] = getgain(tf(50,:),1/200);
    plot(f,gain)
%    set(gca,'xticklabel',num2cell(1./get(gca,'xtick')))
    set(gca,'xticklabel',{'Inf','2000','1000','667','500',400'})
    ylabel('Gain')
    xlabel('Period (years)')
end

% multiply by tv
%test = bsxfun(@times,v230,cSd(:))*(bsxfun(@times,v230,cSdi(:))'*repmat(urtv(:,ii),101,1));

%% Project an interhemispheric pattern
load ../../tmi/tracerobs_4deg_33lev_woce
ip = plotsurf(2*((LAT(jt)>0)' & kt==1)-1);
ip = ip(~isnan(ip))/sqrt(nansum(ip(:).^2)); % normalize to 1

figure(2)
for ii = 1:1
    pii = (ii-1)*3 + 1; % subplot index
    subplot(L,3,pii:(pii+1))
    plotsurf(ps);
    shading flat
    axis equal
    axis tight
    colorbar
    % compute the temporal filter
    ps = ip; % spatial pattern
    
    Us = reshape(ps(:)'*reshape(q*uq,2806,[]),101,230);
    Vs = reshape(ps(:)'*reshape(p*vq,2806,[]),101,230);
    tf = Us*sq*Vs';
  
    % plot the temporal filter gain
    subplot(L,3,pii+2)
    [f,gain] = getgain(tf(50,:),1/200);
    plot(f,gain)
%    set(gca,'xticklabel',num2cell(1./get(gca,'xtick')))
    set(gca,'xticklabel',{'Inf','2000','1000','667','500',400'})
    ylabel('Gain')
    xlabel('Period (years)')
end
test = bsxfun(@times,v230,cSd(:))*(bsxfun(@times,v230,cSdi(:))'*repmat(ip,101,1))

quqr = reshape(q*uq,2806,[]);
pvqr = reshape(p*vq,2806,[]);
%% Now: pointwise
gainmat = [];
for ii = 1:2806
    
    Us = reshape(quqr(ii,:),101,230);
    Vs = reshape(pvqr(ii,:),101,230);
    tf = Us*sq*Vs';
  
    % plot the temporal filter gain
    [f,gain] = getgain(tf(50,:),1/200);
    gainmat(:,ii) = gain;
end

save gainmat f gainmat
