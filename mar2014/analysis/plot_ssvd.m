% plot_ssvd
% run this in the run directory!

%% Setup

clear
close all

k = 230;

addpath ../../data/
addpath ../../../export_fig/
addpath ../utils

TSTEP = 200;
%load g
% if you change stlims, you will have to recompute Gbusv!
stlims = [-25000 -5000];
sig = 0.2; % permil std dev. this is also built into mkD.m

% time axis for the reconstruction
t = stlims(1):TSTEP:stlims(2);
nt = length(t);

ntg = 25; % number of time steps in the gn fn vectors
gam = 1; % tapering parameter

fsp = [1075 1439 1550 821]; % full screen position; for plotting

load ../../tmi/tracerobs_4deg_33lev_woce.mat
% number of surface grid boxes
nr = sum(kt==1);

load ../setup/Dlsfile.mat
%% plot ICs

load icout4
figure(1)
set(1,'color','w')%,'position',[520   378   560   420])
xicn = xic;
xicn(isnan(xicn)) = 0;
plotsurf0(xicn);
axis equal, axis tight
export_fig('-pdf',['Figs/xic_sol',num2str(kic)])

%yic-yold
[u1 s1 v1] = svd(Cxxic);
plotsurf0(u1(:,1));
caxis([-0.03,0.03])
axis equal, axis tight
export_fig('-pdf',['Figs/xic_Cxx',num2str(kic),'eig1'])

plotsurf0(u1(:,2));
axis equal, axis tight
export_fig('-pdf',['Figs/xic_Cxx',num2str(kic),'eig2'])

close
SCALE = 700;
axis equal, axis tight, axis off
% plot the diagonal of the unweighted solution resolution matrix
out = bubbsurf(-diag((diag(ciSd)*(vic(:,1:kic)*vic(:,1:kic)'))*diag(1./ciSd)),10e-4,700);
lr = (max(-out(:)));
scatter(2,-90,SCALE*0.5,'b')
text(15,-90,[num2str(0.5,1) ' permil'])
scatter(71,-90,SCALE*0.1,'b')
text(82,-90,[num2str(0.1,1) ' permil'])
scatter(135,-90,SCALE*0.05,'b')
text(143,-90,[num2str(0.05,1) ' permil'])
scatter(200,-90,SCALE*0.01,'b')
text(207,-90,[num2str(0.01,1) ' permil'])
export_fig('-pdf',['Figs/xic_Tvd',num2str(kic)])

%% plot solution and Cxx

solm = ['solout' num2str(k)];
load([solm])
%Cxxm = ['Cxxout',num2str(k)];
%load([Cxxm])

figure(2)
close
set(gcf,'color','w')
plotys2(Dwd,yestu,Dwdt,core_order,yil,Pm)
export_fig('-pdf',['Figs/yest',num2str(k)])

plot_xest_ts_novb(k)