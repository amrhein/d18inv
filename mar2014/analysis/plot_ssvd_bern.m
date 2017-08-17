% plot_ssvd
% run this in the run directory!

%% Setup

clear
close all

addpath ../../data/
addpath ../utils
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

%% plot solution and Cxx for all cases of k
k = 230;

solm = ['solout' num2str(k)];
load([solm])
Cxxm = ['Cxxout',num2str(k)];
load([Cxxm])

core_order_old = {
    'GeoB1711'
    'GeoB9526_4'
    'M35003_4'
    'MD07_3076'
    'MD98_2165'
    'MD99_2334K'
    'NIOP_905'
    'TR163_31b' };

figure(1)
clf
set(gcf,'color','w')
plotys2(Dwd,yestu,Dwdt,core_order_old,yil,Pm)
export_fig('-eps',['Figs/yest_bern_',num2str(k)])

figure(3)
clf
set(gcf,'color','w')
plotys2(Dwd,yestu*NaN,Dwdt,core_order_old,yil,Pm)
export_fig('-eps',['Figs/yest_bern_only_data',num2str(k)])


figure(2)
clf
load coast
hold on
plot(long,lat,'color',0.7*ones(1,3))
axis equal
axis tight
set(gcf,'color','w')

bdir = '/net/confucius/raid13/home/amrhein/d18inv/data/benthic';
bcores = dir([bdir '/*.mat']);

lbc = length(bcores)
for ii = 1:lbc
    eval(['load ' bdir '/' bcores(ii).name ';']);
    name = strrep(bcores(ii).name,'.mat','');
    eval(['lat=' name '.lat' ';'])
    eval(['lon=' name '.lon' ';'])
    eval(['depth=' name '.depth' ';'])

    plot(lon,lat,'ko','markerfacecolor','k')
    text(lon+2,lat-3,core_order_old(ii),'color','k','fontweight', ...
         'bold','interpreter','none','fontsize',14);
    text(lon+2,lat-10,[num2str(depth), ' m'],'color','k','fontweight', ...
         'bold','interpreter','none','fontsize',14);
end

axis off

set(gcf,'position',[1075         315        1024         482])

export_fig('-eps','Figs/bent_map_bern')