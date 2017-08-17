%function [] = plot_xest_ts_Tv(k)p
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
bind = indmap(~isnan(indmap) & (latmg<-40 & latmg>-60 ) & (lonmg>=30 & lonmg<=140));
reginds(:,2) = 0;
reginds(bind,2) = 1;

bind = indmap(~isnan(indmap) & (latmg<-65 ) & (lonmg>=160 & lonmg<=200));
reginds(bind,3) = 1;

bind = indmap(~isnan(indmap) & (latmg<-52 & latmg>-72 ) & (lonmg>=220 & lonmg<=280));
reginds(bind,4) = 1;

bind = indmap(~isnan(indmap) & (latmg<-60) & (lonmg>=290 & lonmg<=330));
reginds(bind,5) = 1;
sb = indmap(~isnan(indmap) & (latmg<-50));
save sb sb % used in nhavgswd

reg_order = {'1' '2' '3' '4' '5'};
%reg_order = {'1' '2' '3' '4'};

nr = length(reg_order);
%
regmap = zeros(1,2806);%plotsurf(0*var(xestmat)');
for ii = 1:nr
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
cmap(1,:) = 0.5*[1 1 1];
cmap(2,:) = [0.9 0.9 0.9];
colormap(cmap)

text(292,70,'1','fontsize',16,'color','b')
text(23,76,'1','fontsize',16,'color','b')
text(40,-45,'2','fontsize',16,'color','b')
text(191,-68,'3','fontsize',16,'color','b')
text(225,-70,'4','fontsize',16,'color','b')
text(326,-65,'5','fontsize',16,'color','b')
axis equal, axis tight
set(gca,'fontsize',16)
load gainmat
SCALE = 700;

hold on
bubbsurf_nomask_k(-gainmat(1,:),10e-4,SCALE);

scatter(16,-80,SCALE*0.5,'k','filled')
text(28,-78,[num2str(0.5,1)],'color','k')
scatter(54,-79,SCALE*0.1,'k','filled')
text(61,-78,[num2str(0.1,1) ],'color','k')
scatter(85,-79,SCALE*0.05,'k','filled')
text(90,-78,[num2str(0.05,1) ],'color','k')
scatter(118,-79,SCALE*0.01,'k','filled')
text(123,-78,[num2str(0.01,1) ],'color','k')

plot([2 30],[-50 -50],'k--')
plot([140 358],[-50 -50],'k--')

box on
axis on

set(gca,'ytick',[-80,0,80])
set(gca,'xtick',[1,90,180,270,359])

format_ticks(gca,{'0^\circ','90^\circ E','180^\circ E','90^\circ W','0^\circ'},...
    {'80^\circ S','0^\circ','80^\circ N'})

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
        text(lonw-30,lat+3,['  ' name],'interpreter','none','fontsize',12,'color','r','fontweight','bold')
    else
        text(lonw,lat+3,['  ' name],'interpreter','none','fontsize',12,'color','r','fontweight','bold')
    end

    pii = pii+1;
end


export_fig('-pdf',['Figs/xest_map_agu2',num2str(k)])
export_fig('-eps',['Figs/xest_map_agu2',num2str(k)])

%% spag plots in a different style
figure(2)
close
figure(2)
h1 = axes;
RS =4; % amount by which regions are staggered
ofs = round((RS*((numel(reg_order)):-1:0)));

numb = 10; % number of resolution bins
mg = min(gainmat(1,(gainmat(1,:)>0)));
cl = [log10(0.01),log10(max(gainmat(1,:)))];
resbincs = (logspace(cl(1),cl(2),numb));
cmat = gray(round(numb*1.25));
cmat = flipud(cmat(1:numb,:));

for ii = 1:length(reg_order)
    offset = ofs(ii);
    hold on
    res = gainmat(1,reginds(:,ii)); % resolution
    all_nonzero =xestmat(:,reginds(:,ii));
    zinds = ~sum(all_nonzero);
    all_nonzero(:,zinds)=[];
    all_nonzero_res = res;
    all_nonzero_res(zinds)=[];

    [~,sind] = sort(all_nonzero_res);
    
    plot(time,0*time+offset,'k','linewidth',2) % zero line

    for jj = 1:numel(all_nonzero_res)
        [~,bin] = min(abs(resbincs-(all_nonzero_res(sind(jj)))));
        plot(time,all_nonzero(:,sind(jj))+offset,'color',cmat(bin,:),'linewidth',1.5)
    end

    text(-24e3,ofs(ii)+2,[ num2str(ii)],'fontsize',16)
end

colormap(cmat)
h = colorbar('southoutside');

caxis(cl)

set(h,'xtick',linspace(cl(1),cl(2),5))
xticko = get(h,'xtick');
set(h,'xticklabel',num2str(10.^xticko',1))
set(get(h,'xlabel'),'String', 'Pointwise resolution','fontsize',16);

%set(gcf,'color','w','position',[ 1076         238         560         559])

set(gcf,'color','w','position',[ 1184          65         416         719])
YSPACE = 1;
gridxy(min(xlim):1000:max(xlim),min(ylim):YSPACE:max(ylim),'linestyle',':','color',[.5 .5 .5])
yl = ylim;
set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})
xlabel('Time (years)','fontsize',16)
ylabel('\delta^1^8O_c anomalies (permil PDB)','fontsize',16)
hold on
set(h1,'ytick',sort([ofs,ofs-1,ofs+1]))
set(h1,'ytickLabel',repmat([-1 0 1],1,numel(get(h1,'ytick'))))
%set(h1,'ytick',fliplr(ofs))
%set(h1,'ytickLabel',zeros(1,numel(get(h1,'ytick'))))
%ylim([31 93])
set(gca,'fontsize',16)
fn = ['Figs/xest_spag_agu',num2str(k)];
export_fig('-eps',fn)

%%



