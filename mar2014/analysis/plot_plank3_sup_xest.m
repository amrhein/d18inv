% a script to plot the planktonic data

clear
cores = dir('../../data/plank/*.mat');
load solout230
xestmat = reshape(xest,2806,[]);
load ../../tmi/tracerobs_4deg_33lev_woce

LIMS = [-25000,-5000];

figure(1)
clf
hold on

figure(3)
clf
hold on

figure(2)
clf
load coast
hold on
plot(long,lat,'g')
axis equal
axis tight
set(gcf,'color','w')
piiNH = 1; % subplot counter
piiSH = 1; % subplot counter

SP = 3; % spacing between records
dmv = [];

% rank by latitude
latv = [];
for ii = 1:length(cores)
    eval(['load ../../data/plank/' cores((ii)).name ';']);
    name = strrep(cores((ii)).name,'.mat','')
    eval(['latv(ii)=' name '.lat' ';'])
end
[~,s_ind] = sort(latv,'descend');

% initialize interp'd time series matrix
tp = -25000:200:-5000;
ltp = length(tp);
lc = length(cores);

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
%    dm = nanmean(data(:,5));
    dm = nanmean(data([1 end],5));
    
    
        figure(1)
        pii = piiNH;
        piiNH = piiNH+1;
        dmvNH(ii) = dm;
        hold all
    plot(data(:,2),data(:,5)-dm-pii*SP,'k')
    plot(-25000:200:-5000,xestmat(tmi_ind,:)-pii*SP-mean(xestmat(tmi_ind,:)),'b')
    text(-24.5e3,-pii*SP+1,name,'interpreter','none','fontsize',10)
    xlim(LIMS)
    %    legend('raw','interp')
    
    figure(2)
    plot(lon,lat,'r*')
    text(lon,lat,['  ' name],'interpreter','none')
    
    pii = pii+1;
end

addpath ../../../export_fig/
figure(1)
set(gcf,'color','w','position',[1774         272         393 502])
set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})
xlabel('time (yrs)')
ylabel('\delta^1^8O_c (permil)')
%export_fig('-pdf','Figs/planks')
mvec = flipud(dmv(:));
%set(gca,'yticklabel',{'0','0','0','0','0','0','0','0','0'});
set(gca,'yticklabel',{'';num2str(mvec,3);''});
ylabel('\delta^1^8O_c')
xlabel('Time (years)')
yl = ylim;
set(gca,'ytick',yl(1):1:yl(2),'yticklabel','');
grid
export_fig('-pdf','Figs/plankall_x')


