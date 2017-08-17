function [] = plot_xest_ts_novb(k)
% plot_xest_ts
% a script to plot time series from solutions to the IBC problem
% no vb = no variance bubble plto
%clear

addpath ../../../export_fig

%   k = 230;
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
%reginds = false(2806,5);
reginds = false(2806,4);

% Define the Northern box
nb = indmap(~isnan(indmap) & (latmg>=46 & latmg<80) & ((lonmg>=0 & lonmg<=30) |...
      (lonmg>=290 & lonmg<=360)));
 reginds(nb,1 ) = 1;

% WPAC
ross = indmap(~isnan(indmap) & latmg>=-10 & latmg<=10 & lonmg>=100 & lonmg<=150);
reginds(ross,2) = 1;

% EPAC
ross = indmap(~isnan(indmap) & latmg>=-10 & latmg<=10 & lonmg>=210 & lonmg<=260);
reginds(ross,3) = 1;

% Define the Southern box
%sb = indmap(~isnan(indmap) & (latmg>=-84 & latmg<-46) & ((lonmg>=0 & lonmg<=30) |...
%      (lonmg>=290 & lonmg<=360)));
sb = indmap(~isnan(indmap) & latmg<-60);
%sb = indmap(~isnan(indmap) & (latmg>=-84 & latmg<-66) & ((lonmg>=0 & lonmg<=30) |...
%      (lonmg>=290 & lonmg<=360)));
 reginds(sb,4) = 1;

% Ross -70:-82 lat, 166:194 lon
%sb = indmap(~isnan(indmap) & latmg>=-82 & latmg<=-70 & lonmg>=166 & lonmg<=194);
%reginds(sb,5) = 1;

%reg_order = {'1' '2' '3' '4' '5'};
reg_order = {'1' '2' '3' '4'};

nr = length(reg_order);
%%
regmap = zeros(1,2806);%plotsurf(0*var(xestmat)');
for ii = 1:nr
    regmap(reginds(:,ii)) = ii+1;
end

close all
rmp = plotsurf(regmap);
set(gcf,'color','w')
pcolor(lonmg,latmg,rmp)
shading flat
cmap = lines(nr+5);
cmap(end,:) = [];
cmap(1,:) = [0.8 0.8 0.8];
colormap(cmap)

text(316,57,'1','fontsize',16,'color','w')
text(112,6,'2','fontsize',16,'color','w')
text(220,6,'3','fontsize',16,'color','w')
%text(170,-70,'5','fontsize',16,'color','w')
text(180,-72,'4','fontsize',16,'color','w')
axis equal, axis tight
export_fig('-pdf',['Figs/xest_map',num2str(k)])

%% spag plots in a different style
figure(1)
clf
h1 = axes;
RS =4; % amount by which regions are staggered
ofs = round((RS*((numel(reg_order)):-1:0)));
for ii = 1:length(reg_order)
    offset = ofs(ii);
    hold on
    plot(time,xestmat(:,reginds(:,ii))'+offset)
    text(-24e3,ofs(ii)+2,[ num2str(ii)],'fontsize',16)
end

%set(gcf,'color','w','position',[ 1076         238         560         559])
set(gcf,'color','w','position',[ 1184          65         416         719])
YSPACE = 1;
gridxy(min(xlim):1000:max(xlim),min(ylim):YSPACE:max(ylim),'linestyle',':','color',[.5 .5 .5])
set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})
xlabel('Time (years)')
ylabel('Permil \delta')
hold on
set(h1,'ytick',fliplr(ofs))
set(h1,'ytickLabel',zeros(1,numel(get(h1,'ytick'))))
%ylim([31 93])
fn = ['Figs/xest_spag',num2str(k)];
export_fig('-pdf',fn)

%%



