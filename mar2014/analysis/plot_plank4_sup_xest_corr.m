% a script to plot the planktonic data

clear
close all
cores = dir('../../data/plank/*.mat');
load solout230
xestmat = reshape(xest,2806,[]);
load ../../tmi/tracerobs_4deg_33lev_woce

LIMS = [-25000,-5000];

figure(1)
clf
hold on
set(gcf,'color','w','position',[1 40 390 761])

figure(3)
clf
hold on
set(gcf,'color','w','position',[1 40 390 761])
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

load gainmat

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

    % % This part interpolates the record to 200 year spacing in
    % % order to compute correlations
 
  % all times and data
   ta = data(:,2);
   ya = data(:,5);

    % eliminate multiple measurements at the same time (the result of
    % samples at the same depth) by using the average
    
    allreps = [];
    for i2 = 1:length(ta)
        repinds = find(ta==ta(i2));
        ya(repinds) = mean(ya(repinds));
        if numel(repinds)>1, allreps = [allreps;repinds(2:end)]; end
    end
    ya(unique(allreps)) = [];
    ta(unique(allreps)) = [];

   tdgi = (ta>=-26000 & ta<-4000);
   t = ta(tdgi);
   y = ya(tdgi);

   tp = -26000:-4000;
   yi = interp1(t,y,tp);

   yis = conv(yi,ones(200,1)/200,'same'); %smoothed
   yiss = yis(1:200:end); % subsampled

   yii = yiss(6:(end-5));

   % exception for MD95-2010, which ends early but is a great
   % record - to compute the mean, pretend that the last value
   % persists until the end of the reconstruction interval
   if strcmp(name,'MD95_2010')
       yiim = yii-(nansum(yii)+sum(isnan(yii))*yii(find(isnan(yii),1)-1))/length(yii);
   else
       yiim = yii-nanmean(yii);
   end
   % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
   
   % Actually compute the correlation
   if ~any(xestmat(tmi_ind,~isnan(yiim)))
      rp = 0;
   else
       rpm = corrcoef(xestmat(tmi_ind,~isnan(yiim)),yiim(~isnan(yiim)));
       rp = rpm(1,2);
   end
   disp(rp)
   if isnan(rp)
       keyboard
   end
   % Actually compute the correlation

%    dm = nanmean(data(:,5));
    dm = nanmean(data([1 end],5));
    
    if lat>0
        figure(1)
        pii = piiNH;
        piiNH = piiNH+1;
        dmvNH(ii) = dm;
    else
        figure(3)
        pii = piiSH;
        piiSH = piiSH+1;
        dmvSH(ii) = dm;
    end
    subplot(1,4,1:3)
    hold on
        plot(data(:,2),data(:,5)-dm-pii*SP,'k')
    %    Checking that the correlation calculations are on sane quantities:
    %    plot(-25000:200:-5000,yiim-pii*SP,'k')
    plot(-25000:200:-5000,xestmat(tmi_ind,:)-pii*SP-mean(xestmat(tmi_ind,:)),'b')
            text(-24.5e3,-pii*SP-0.5,[name],'interpreter','none','fontsize',10)
        text(-4.5e3,-pii*SP-0.5,['$r$= ',num2str(gainmat(1,tmi_ind),1)],'interpreter','latex','fontsize',10)
        text(-4.5e3,-pii*SP+0.5,['$R$ = ',num2str(rp,2)],'interpreter','latex','fontsize',10)
    % indicate gain
    xlim(LIMS)
    %    legend('raw','interp')

    figure(2)
    plot(lon,lat,'r*')
    text(lon,lat,['  ' name],'interpreter','none')
    
    pii = pii+1;
end

%%

addpath ../../../export_fig/
figure(1)
%set(gcf,'color','w','position', [1774         272         393 502])
set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})
xlabel('Time (years)')
ylabel('\delta^1^8O_c anomalies (permil PDB)')
%export_fig('-pdf','Figs/planks')
mvec = flipud(dmvNH(:));
%set(gca,'yticklabel',{'';num2str(mvec,3);''});
%set(gca,'yticklabel',{num2str(mvec,3)});
set(gca,'ytick',SP*(-length(dmvNH):-1))
set(gca,'yticklabel',{'0','0','0','0','0','0','0'});
ylabel('\delta^1^8O_c anomalies (permil VPDB)')
xlabel('Time (years)')
%yl = ylim;
%set(gca,'ytick',yl(1):1:yl(2),'yticklabel','');
ylim([-24 -1])
gridxy(min(xlim):2500:max(xlim),min(ylim):max(ylim),'linestyle',':','color',[.5 .5 .5])
export_fig('-pdf','Figs/plankNH_x_corr')

figure(3)
%set(gcf,'color','w','position', [1774         272         393 502])
set(gcf,'color','w','position',[1 40 390 761])
set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})
xlabel('Time (years)')
set(gca,'ytick',SP*(-sum(~~dmvSH):-1))
set(gca,'yticklabel',{'0','0','0','0','0','0'});
ylabel('\delta^1^8O_c anomaly (permil VPDB)')
gridxy(min(xlim):2500:max(xlim),min(ylim):max(ylim),'linestyle',':','color',[.5 .5 .5])
ylim([-20 -1])

export_fig('-pdf','Figs/plankSH_x_corr')
