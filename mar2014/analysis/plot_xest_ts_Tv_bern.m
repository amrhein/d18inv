%function [] = plot_xest_ts_Tv(k)

clear
close all
k = 230

addpath ../../../export_fig

n = 25;
time = -25000:200:-5000;

load ../../tmi/tracerobs_4deg_33lev_woce.mat

figure(1)
[lonmg,latmg] = meshgrid(2:4:358,-90:4:88);

load gainmat
SCALE = 700;

hold on
bubbsurf_k(-gainmat(1,:),10e-4,SCALE);
axis off
scatter(2,-90,SCALE*0.5,'k')
text(15,-90,[num2str(0.5,1)],'color','k')
scatter(71,-90,SCALE*0.1,'k')
text(82,-90,[num2str(0.1,1) ],'color','k')
scatter(135,-90,SCALE*0.05,'k')
text(143,-90,[num2str(0.05,1) ],'color','k')
scatter(200,-90,SCALE*0.01,'k')
text(207,-90,[num2str(0.01,1) ],'color','k')


export_fig('-eps',['Figs/tv_map_bern',num2str(k)])


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
    %    if lonw>180
    %    text(lonw-30,lat+3,['  ' name],'interpreter','none','fontsize',14,'color','r','fontweight','bold')
    %else
    %    text(lonw,lat+3,['  ' name],'interpreter','none','fontsize',14,'color','r','fontweight','bold')
    %end

    pii = pii+1;
end

%set(gcf,'position',[1075         165        1236         632])
set(gca,'fontsize',16)

axis equal, axis tight
box on
axis on

set(gca,'ytick',[-80,0,80])
set(gca,'xtick',[1,90,180,270,359])

format_ticks(gca,{'0^\circ','90^\circ E','180^\circ E','90^\circ W','0^\circ'},...
    {'80^\circ S','0^\circ','80^\circ N'})


export_fig('-eps',['Figs/tv_map_bern_plank',num2str(k)])

