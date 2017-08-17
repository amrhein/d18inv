%function [] = plot_xest_ts(xest,n,region)
% plot_xest_ts
% a script to plot time series from solutions to the IBC problem
% 19 Feb 2014
% D Amrhein
clear

addpath ../../../export_fig
load ../../tmi/Vtot_DEA_30Sept.mat
v=(Vtot(~isnan(Vtot)));

load solout230.mat
n = 25;
region = 1; %globe

load ../../tmi/c_all_4deg.mat
load ../../tmi/tracerobs_4deg_33lev_woce.mat

% N = 25;
xestmat = reshape(xest,2806,[])';
% subselect by region
regs = logical(c_all(kt==1,region));
[val,ind] = sort((var(xestmat)),'descend');
is_top = var(xestmat)'>=val(n);

figure(2)
clf
bubbsurf_nomask(c_all(kt==1,region).*var(xestmat)',0.001,100);

[lonmg,latmg] = meshgrid(2:4:358,-90:4:88);

%% Define some regions

indmap = plotsurf(0*var(xestmat)');
indmap(~isnan(indmap)) = 1:2806;
reginds = false(2806,7);

% Gnland: 50:62 lat, 306:326 lon
gnl = indmap(~isnan(indmap) & latmg>=54 & latmg<64 & lonmg>=306 & lonmg<=346);
reginds(gnl,1) = 1;

% GIN1: 58:78 lat, 334:360 lon and 0:10 lon 
gin1 = indmap(~isnan(indmap) & (latmg>=64 & latmg<=78 & ...
    lonmg>=342 & lonmg<=360) |...
    (latmg>=68 & latmg<=78 & lonmg<=2));
gin1 = [gin1;2768;2800];
reginds(gin1,2) = 1;

% Med: 18:34 lat, 10:38 lon
med = indmap(~isnan(indmap) & (latmg>=18 & latmg<54) & ((lonmg>=0 & lonmg<=38) |...
     (lonmg>=346 & lonmg<=360)));
reginds(med,3) = 1;

% GIN2 (east)
gin2 = indmap(~isnan(indmap) & (latmg>=60 & latmg<=78 & ...
    (lonmg>=0 & lonmg<=20)) &...
    ~(latmg>=68 & latmg<=78 & lonmg<=2));
reginds(gin2,4) = 1;

% GIN3 (south)
% gin3 = indmap(~isnan(indmap) & latmg>=54 & latmg<=64 & ...
%     ((lonmg>=350 & lonmg<=360)));
% gin3(gin3==2800 | gin3 == 2768)=[];
% 
% reginds(gin3,7) = 1;

% ARC: 70:82 lat, 34:62 lon
arc = indmap(~isnan(indmap) & latmg>=70 & latmg<82 & lonmg>=34 & lonmg<=62);
reginds(arc,5) = 1;

% Aus: -38:-54 lat, 66:138 lon
aus = indmap(~isnan(indmap) & latmg>=-54 & latmg<=-38 & lonmg>=66 & lonmg<=138);
reginds(aus,6) = 1;

% Ross -70:-82 lat, 166:194 lon
ross = indmap(~isnan(indmap) & latmg>=-82 & latmg<=-70 & lonmg>=166 & lonmg<=194);
reginds(ross,7) = 1;

% Weddell: -58:-78 lat, 298:334 lon
%wed = indmap(~isnan(indmap) & latmg>=-78 & latmg<=-78 & lonmg>310 & lonmg<=334);
%reginds(wed,3) = 1;

% SO2: -70:-52 lat, 222:290 lon
so2 = indmap(~isnan(indmap) & (latmg>=-78 & latmg<-52) & ((lonmg>=222 & lonmg<=326)));
reginds(so2,8) = 1;

%reg_order = {'aus' 'ross' 'so2' 'gnl' 'gin1' 'gin2' 'arc' 'med'};
reg_order = {'1' '2' '3' '4' '5' '6' '7' '8'};
nr = length(reg_order);
%%
regmap = zeros(1,2806);%plotsurf(0*var(xestmat)');
for ii = 1:nr
    regmap(reginds(:,ii)) = ii+1;
end

figure(2)
rmp = plotsurf(regmap);
clf
set(gcf,'color','w')
pcolor(lonmg,latmg,rmp)
shading flat
cmap = lines(nr+5);
cmap(end,:) = [];
cmap(1,:) = [0.8 0.8 0.8];
colormap(cmap)

hold on
%plotsurf(c_all(kt==1,region).*var(xestmat')');
SCALE = 50;
bubbsurf_nomask(-var(xestmat)',0.001,SCALE); % minus is for color
set(gca,'ydir','n')

text(306,57,'1','fontsize',16,'color','w')
text(342,78,'2','fontsize',16,'color','w')
text(342,46,'3','fontsize',16,'color','w')
text(8,40,'3','fontsize',16,'color','w')
text(18,78,'4','fontsize',16,'color','w')
text(34,82,'5','fontsize',16,'color','w')
text(66,-38,'6','fontsize',16,'color','w')
text(166,-70,'7','fontsize',16,'color','w')
text(222,-62,'8','fontsize',16,'color','w')


lr = (max(abs(xestmat(:))));
scatter(2,-90,SCALE*lr,'b')
text(15,-90,[num2str(lr,1) ' permil'])
scatter(71,-90,SCALE*5,'b')
text(82,-90,[num2str(5,1) ' permil'])
scatter(135,-90,SCALE*1,'b')
text(143,-90,[num2str(1,1) ' permil'])
scatter(195,-90,SCALE*0.1,'b')
text(202,-90,[num2str(0.1,1) ' permil'])
axis off

export_fig('-pdf','Figs/xest_map')

%set(gcf,'position',[154   598   918   512])

%%
figure(1)
clf
time = -25000:200:-5000;
for ii = 1:length(reg_order)
    subplot(2,4,ii)
    %eval(['plot(time,xestmat(:,' reg_order{ii} '))'])
    plot(time,xestmat(:,reginds(:,ii))')
    title(reg_order(ii))
    vr = v(reginds(:,ii));
    vn = vr/sum(vr);
%    hold on
%    plot(time,xestmat(:,reginds(:,ii))*vn,'k','linewidth',2)
    axis tight
    ylim([-10 10])
end
%% spag plots in a different style
figure(1)
clf
time = -25000:200:-5000;
h1 = axes;
RS = 11; % amount by which regions are staggered
ofs = 2*round((RS*((numel(reg_order)):-1:0)+...
    flipud(cumsum([0 2 5 -1 7 5 3 1 0])))/2);
for ii = 1:length(reg_order)
    %offset = RS*(numel(reg_order)+1-ii)+ofs(ii)
    offset = ofs(ii);
    hold on
%    if ii>1
%        offset = max(max(xestmat(:,reginds(:,ii)))-min(min(xestmat(:,reginds(:,ii-1))));
%    end
    plot(time,xestmat(:,reginds(:,ii))'+offset)
    %text(-24e3,RS*ii+4,[ num2str(numel(reg_order)+1-ii)],'fontsize',16)
    text(-24e3,ofs(ii)+2,[ num2str(ii)],'fontsize',16)
end
set(gcf,'color','w','position',[ 1076         238         560         559])
YSPACE = 2;
gridxy(min(xlim):1000:max(xlim),min(ylim):YSPACE:max(ylim),'linestyle',':','color',[.5 .5 .5])
set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})
xlabel('Time (years)')
ylabel('Permil \delta')
hold on
% h2 = axes('position',get(h1,'position'),...
%     'xtick',get(h1,'xtick'),...
%     'xticklabel',get(h1,'xticklabel'),...
%     'yaxislocation','right','color','none');
%set(h1,'ytick',min(ylim(h1)):RS:max(ylim(h1)))
set(h1,'ytick',fliplr(ofs))
set(h1,'ytickLabel',zeros(1,numel(get(h1,'ytick'))))
ylim([31 93])
export_fig('-pdf','Figs/xest_spag')

%%

figure(3)
bubbsurf(-v/sum(v),0,5e3);
set(gcf,'position',[154   598   918   512])


%% 
figure(4)
vx = var(xestmat);
v=(Vtot(~isnan(Vtot)));
scatter((v/sum(v)),(vx),'.k')
set(gca,'xscale','log','yscale','log')
%ylim([1e-16,1e2])

xlabel('Ocean-filling fraction')
ylabel('Variance of inverse solution')
set(gcf,'color','w')

% %%
% figure(5)
% clf
% load /home/dan/Dropbox/2013-2014/Amrhein_deglacial_d18O/fweddell/resout413.mat
% tt = diag(Tv(2*2806+1:3*2806,:));
% %vx = var(xestmat);
% v=(Vtot(~isnan(Vtot)));
% %scatter((v/sum(v)),vx,'.k')
% scatter((v/sum(v)),tt,'.k')
% %scatter(vx,tt,'.k')
% set(gca,'xscale','log','yscale','log')
% ylim([1e-16,1])
% set(gca,'ytick',[1e-16,1e-12,1e-8,1e-4,1])
% 
% xlabel('Ocean-filling fraction')
% ylabel('Resolution')
% set(gcf,'color','w','position',[1850         228         560         249])
% export_fig('-pdf','Figs/res-fill-scatter')
