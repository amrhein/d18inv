%function [] = plot_xest_ts_Tv(k)
% plot_xest_ts
% a script to plot time series from solutions to the IBC problem
% no vb = no variance bubble plto
%clear

clear
k = 230

addpath ../../../export_fig

load(['solout',num2str(k)])
n = 25;
region = 1; %globe
time = -25000:200:-5000;

load ../../tmi/tracerobs_4deg_33lev_woce.mat

% N = 25;
xestmat = reshape(xest,2806,[])';


figure(2)
[lonmg,latmg] = meshgrid(2:4:358,-90:4:88);

%% Define some regions

indmap = plotsurf(0*var(xestmat)');
indmap(~isnan(indmap)) = 1:2806;
reginds = false(2806,5);

% Define the Northern box
%nb = indmap(~isnan(indmap) & (latmg>=46 & latmg<80) & ((lonmg>=0 & lonmg<=30) |...
%      (lonmg>=300 & lonmg<=360)));

nb = indmap(~isnan(indmap) & (latmg>=46 & latmg<80) & ((lonmg>=0 & ...
                                                  lonmg<=30) | (lonmg>=290 & lonmg<=360)));
%nb = indmap(~isnan(indmap) & (latmg>=50 & latmg<80) & (lonmg>=300 & lonmg<=360));
reginds(nb,1 ) = 1;

% Southern boxes: uneven
bind = indmap(~isnan(indmap) & (latmg<-50));
reginds(:,2) = 0;
reginds(bind,2) = 1;

reg_order = {'1' '2' '3' '4' '5'};
%reg_order = {'1' '2' '3' '4'};

nr = length(reg_order);
%
regmap = zeros(1,2806);%plotsurf(0*var(xestmat)');
for ii = 1:2
    %        regmap(reginds(:,ii)) = ii+1;
        regmap(reginds(:,ii)) = 2;
end

close all

rmp = plotsurf(regmap);
set(gcf,'color','w')
pcolor(lonmg,latmg,rmp)
shading flat
%cmap = lines(nr+5);
%cmap(end,:) = [];
cmap = [];
cmap(1,:) = [0.8 0.8 0.8];
cmap(2,:) = [0.8 0.8 0.8];
colormap(cmap)

%text(292,70,'NA','fontsize',16,'color','r')
%text(23,76,'NA','fontsize',16,'color','r')
%text(40,-45,'2','fontsize',16,'color','r')
%text(191,-68,'3','fontsize',16,'color','r')
%text(225,-70,'4','fontsize',16,'color','r')
%text(326,-65,'SO','fontsize',16,'color','r')
axis equal, axis tight
load gainmat
SCALE = 700;

hold on
%bubbsurf_nomask_k(-gainmat(1,:),10e-4,SCALE);


scatter(2,-90,SCALE*0.5,'k')
text(15,-90,[num2str(0.5,1)],'color','k')
scatter(71,-90,SCALE*0.1,'k')
text(82,-90,[num2str(0.1,1) ],'color','k')
scatter(135,-90,SCALE*0.05,'k')
text(143,-90,[num2str(0.05,1) ],'color','k')
scatter(200,-90,SCALE*0.01,'k')
text(207,-90,[num2str(0.01,1) ],'color','k')

%{
 lr = 1;
 scatter(2,-90,SCALE*1,'b')
 text(15,-90,[num2str(lr,1) ])
 scatter(71,-90,SCALE*lr/2,'b')
 text(82,-90,[num2str(lr/2,1)])
 scatter(135,-90,SCALE*0.1,'b')
 text(143,-90,[num2str(0.1,1) ])
 scatter(195,-90,SCALE*0.01,'b')
 text(202,-90,[num2str(0.01,1)])
 %}

 %plot([2 30],[-50 -50],'k--')
 %plot([140 358],[-50 -50],'k--')

% from plot_plank_3

cores = dir('../../data/plank/*.mat');
lc = length(cores)

hold on

piiNH = 1; % subplot counter
piiSH = 1; % subplot counter

dmv = [];

% rank by latitude
latv = [];
for ii = 1:lc
    eval(['load ../../data/plank/' cores((ii)).name ';']);
    name = strrep(cores((ii)).name,'.mat','')
    eval(['latv(ii)=' name '.lat' ';'])
end
[~,s_ind] = sort(latv,'descend');
pii = 0;

for ii = 1:lc
    % load in order by sorted lat
    eval(['load ../../data/plank/' cores(s_ind(ii)).name ';']);
    name = strrep(cores(s_ind(ii)).name,'.mat','')
    eval(['species=' name '.species' ';'])
    eval(['data=' name '.data' ';'])
    eval(['lat=' name '.lat' ';'])
    eval(['lon=' name '.lon' ';'])
    eval(['depth=' name '.depth' ';'])
    
    if lon<0
        lonw = lon+360; 
    else
        lonw=lon;
    end
    
    % surface only
    jts = jt(kt==1);
    its = it(kt==1);
    LATs = LAT(jts);
    LONs = LON(its);
    
    [~,tla] = min(abs(LATs-lat));
    [~,tlo] = min(abs(LONs-lonw)); % wraparound issues?
    tmi_ind = find(( LONs==LONs(tlo) & LATs==LATs(tla)),1)

    % if this is not a valid location, find the nearest long
    searchind = 2;
    while isempty(tmi_ind)
        [yy,iy] = sort(abs(LONs-lonw)); % wraparound issues?
        tmi_ind = find(( LONs==LONs(iy(searchind)) & LATs==LATs(tla)),1)
        searchind= searchind+1;
    end

    data(:,2) = -data(:,2);
    data = flipud(data);
    if all(isnan(data(:,2)))
        continue
    elseif max(abs(data(:,2)))<1000
        data(:,2) = 1000*data(:,2);
    end
    
    data(~~sum(data(:,5)<-8,2),:) = nan;
    dm = nanmean(data([1 end],5));
    
    plot(lonw,lat,'rsq','markerfacecolor','r')
        if lonw>180
        text(lonw-30,lat+3,['  ' name],'interpreter','none','fontsize',14,'color','r','fontweight','bold')
    else
        text(lonw,lat+3,['  ' name],'interpreter','none','fontsize',14,'color','r','fontweight','bold')
    end

    pii = pii+1;
end

%set(gcf,'position',[1075         165        1236         632])
axis equal, axis tight

axis off
export_fig('-eps',['Figs/xest_map_bern5',num2str(k)])
