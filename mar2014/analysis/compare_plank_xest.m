% a script to plot the planktonic data and compare to xest from the benthic
% inversion

clear
close all

xestmatr = sum(xestmatv(regs,:));

cores = dir('../Data/plankd18o/*.mat');

LIMS = [-25000,-5000];

figure(1)

figure(2)
load coast
hold on
plot(long,lat,'g')
    axis equal
    axis tight
set(gcf,'color','w')
title('Locations of planktonic cores')

for ii = 1:length(cores)
    % get species info
    eval(['load ../Data/plankd18o/' cores(ii).name ';']);
    name = strrep(cores(ii).name,'.mat','');
    eval(['species=' name '.species' ';'])
    eval(['data=' name '.data' ';'])
    eval(['lat=' name '.lat' ';'])
    eval(['lon=' name '.lon' ';'])
    eval(['depth=' name '.depth' ';'])
    
    data(:,2) = -data(:,2);
    data = flipud(data);
    if max(abs(data(:,2)))<1000
        data(:,2) = 1000*data(:,2);
    end
    
    figure(1)
    subplot(3,4,ii)
data(~~sum(data(:,5)<-8,2),:) = nan;
plot(data(:,2),data(:,5)-nanmean(data(:,5)))
hold all

plot(-25000:200:-5000,(xestmatr)')
    
    xlim(LIMS)
%    legend('raw','interp')
    set(gcf,'color','w','position',[1075         159        1550         821])
    set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})
    title([name,': lat=',num2str(lat),' lon=',num2str(lon),' depth=',num2str(depth), ],'interpreter','none')
    xlabel('time (yrs)')
    ylabel('\delta^1^8O_c (permil)')
    legend('plank','benthic recon')
    
    figure(2)
    plot(lon,lat,'r*')
    if ii == 1 | ii==2
        text(lon,lat,'HU73_031_7','interpreter','none')
    else
    text(lon,lat,['  ' name],'interpreter','none')
    end
    
end

addpath ../../../export_fig/

%figure(1)
%export_fig('-pdf','Figs/planks')
%figure(2)
%set(gcf,'fontsize',10)
export_fig('-pdf','Figs/plank_locs')

