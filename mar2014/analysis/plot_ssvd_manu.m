% plot_ssvd
% run this in the run directory!

%% Setup

clear
close all

addpath ../../data/
addpath ../../../export_fig/

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

load icout3
figure(1)
close
figure(1)
set(1,'color','w')%,'position',[520   378   560   420])
%plotsurf(xic);
SCALE = 50;
bubbsurf(-xic,0.001,SCALE);
axis equal, axis tight
lr = (max(xic(:)));
scatter(2,-90,SCALE*lr,'b')
text(15,-90,[num2str(lr,1) ' permil'])
scatter(71,-90,SCALE*2,'b')
text(82,-90,[num2str(2,1) ' permil'])
scatter(135,-90,SCALE*1,'b')
text(143,-90,[num2str(1,1) ' permil'])
scatter(195,-90,SCALE*0.1,'b')
text(202,-90,[num2str(0.1,1) ' permil'])
axis off
export_fig('-pdf',['Figs/xic_sol',num2str(kic)])


yic-yold

%% plot solution and Cxx for all cases of k
k = 230;

solm = ['solout' num2str(k)];
load([solm])
Cxxm = ['Cxxout',num2str(k)];
load([Cxxm])

figure(1)
clf
set(gcf,'color','w')
plotys2(Dwd,yestu,Dwdt,core_order,yil,Pm)
export_fig('-pdf',['Figs/yest',num2str(k)])

%%
Cxxr = reshape(sqrt(Cxx),2806,[]);
%Cxxrt = Cxxr(1:2806,:);
[uxrt sxrt vxrt] = svd(Cxxr);
%
figure(4)
clf
%set(gcf,'position',[ 1262         266         569         531],'color','w')
SCALE = 500;
bubbsurf(uxrt(:,1)*sxrt(1,1)/10,0.001,SCALE)
axis off
% make a legend
hold on
% largest value, rounded
lr = (max(abs((uxrt(:,1)*sxrt(1,1)/10))));
scatter(2,-90,SCALE*lr,'b')
text(15,-90,[num2str(lr,1) ' permil'])
scatter(71,-90,SCALE*lr/2,'b')
text(82,-90,[num2str(lr/2,1) ' permil'])
scatter(135,-90,SCALE*0.1,'b')
text(143,-90,[num2str(0.1,1) ' permil'])
scatter(195,-90,SCALE*0.02,'b')
text(202,-90,[num2str(0.02,1) ' permil'])
axis off
xlabel(' ')
export_fig('-pdf',['Figs/CxxEOF' num2str(k)])%,'renderer','painters')

figure(5)
clf
plot(t,-vxrt(:,1)*10)%*sxrt(1,1))
svar = diag(sxrt).^2/sum(diag(sxrt.^2));
ylim([-0.2 1.5])
%set(gcf,'position',[53     2   796   236])
set(gca,'xtick',[-25000:2500:-5000])
set(gca,'xticklabel',num2str(get(gca,'xtick')',5))
grid
xlabel('Time (years)')
ylabel('Normalized units')
set(gcf,'position',[1726         522         560         176],'color','w')

export_fig('-pdf',['Figs/CxxPC' num2str(k)])%,'renderer','painters')
%print('-dpdf',['Figs/Cxx' num2str(k)])
%% load res mats

% load fweddell/resout413.mat
% Tvs = Tv((2*nr+1):3*nr,:);
% load /home/dan/Desktop/TMI_v6-1/Vtot_DEA_30Sept.mat
% Vni = Vtot(~isnan(Vtot));
% sv = sum(Vni);

%% Plot resolution matrices

% plot Tv
% 1: variance reduction plot. The variable Tvs has columns corresponding
% to the various surface locations at the time index value of 50. The rows
% of Tvs correspond to all spaces at the times 48, 49, 50, 51.

% figure(1), close, figure(1)
% SCALE = 100;
% bubbsurf(-diag(Tvs),0.0001,SCALE);
% axis equal, axis tight
% set(gcf,'color','w')
% axis off
% 
% lr = 1;
% scatter(2,-90,SCALE*1,'b')
% text(15,-90,[num2str(lr,1) ])
% scatter(71,-90,SCALE*lr/2,'b')
% text(82,-90,[num2str(lr/2,1)])
% scatter(135,-90,SCALE*0.1,'b')
% text(143,-90,[num2str(0.1,1) ])
% scatter(195,-90,SCALE*0.01,'b')
% text(202,-90,[num2str(0.01,1)])
% axis off
% 
% 
% export_fig('-pdf',['Figs/Tv_diag' num2str(413)],1)
% 

%%
% load fweddell/resout413.mat
% close all
% set(gcf,'color','w')
% %tvind = 28;
% %tvind = 30;
% tvind = 1400; % cent pac
% fi = find(kt==1,tvind);
% ltvind = fi(end);
% tvlat = LAT(jt(ltvind));
% tvlon = LON(it(ltvind));
% 
% %cax = [-0.03    0.3];
% 
% addpath ../../export_fig/cm_and_cb_utilities/
% set(gcf,'position',[1009          19         588         778])
% 
% subplot(3,1,1)
% plotsurf(Tv(2806*1+1:2806*2,tvind));
% hold on
% plot(tvlon,tvlat,'go','markersize',8,'linewidth',3)
% colormap(b2r(caxis));
% %caxis(cax/10)
% colorbar
% axis equal, axis tight
% title('T = -15200')
% 
% subplot(3,1,2)
% plotsurf(Tv(2806*2+1:2806*3,tvind));
% hold on
% plot(tvlon,tvlat,'go','markersize',8,'linewidth',3)
% colormap(b2r(caxis));
% colorbar
% axis equal, axis tight
% title('T = -15000')
% 
% subplot(3,1,3)
% plotsurf(Tv(2806*3+1:2806*4,tvind));
% hold on
% plot(tvlon,tvlat,'go','markersize',8,'linewidth',3)
% %caxis(cax/100)
% hc = colorbar;
% axis equal, axis tight
% colormap(b2r(caxis,37));
% title('T = -14800')
% 
% 
% %axis equal, axis tight
% %export_fig('-pdf',['Figs/Tv_lobes' num2str(413)],1)
% 
% %% plot Tu
% %plotys_Tu(Dwd,yestu,Dwdt,core_order,yil,Pm,Tu)
% plotys_Tu2(Dwd,yestu,Dwdt,core_order,yil,Pm,Tu)

%export_fig('-pdf',['Figs/Tu' num2str(413)],1)

%% chi2

% load fweddell/solout300.mat
% var(chol(inv(Pm))*(Dwd-yestu))
% close all
% hist(chol(inv(Pm))*(Dwd-yestu))
% set(gcf,'color','w')
% xlabel('permil \delta^1^8O_c_c (normalized)')
% ylabel('bin counts')

%export_fig('-pdf',['Figs/chi2' num2str(413)],1)


