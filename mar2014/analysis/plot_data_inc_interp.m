clear
close all
addpath ../../export_fig/

LIMS = [-25000,-5000];
TSTEP = 200;
load /home/dan/Dropbox/2012-2013/Data/core_order_new.mat
load /home/dan/Dropbox/2012-2013/Data/tracerobs_4deg_33lev_woce.mat

% get interpolated stuff
load /home/dan/Dropbox/2012-2013/MSfinal/4weddell/Dlsfile.mat
[co] = yv2coreswd(Dwd,yil);
[cop] = yv2coreswd(sqrt(diag(Pm)),yil);
cote = yv2coreswd(Dwdt,yil);
cot = cote;

clf
set(gcf,'color','w')

% for raw data
path1 = '/home/dan/Dropbox/2012-2013/Data/benthicTMI/';
lims = [-25000,-5000];
[D,tD,Di,core_order_new] = mkDni(path1,lims);
files = dir([path1 '*.mat']);
lf = length(files);


for ii = 1:length(core_order)
    % get species info
    eval(['load ' path1 files(ii).name ';']);
    name = strrep(files(ii).name,'.mat','');
    eval(['species=' name '.species' ';'])
    
    % get raw data
    ly = sum(Di==ii);
    y = D(Di==ii);
    t = tD(Di==ii);
    N = 0.2^2*eye(ly);

    ly = length(y);
    ym = get_ym(y,t,N);
    %ym = d18mv;
    ymm = y-ym;

    %subplot(1,3,3)
    hold on
%    plot(t,d18-mean(d18) + ii*2,'.')

SP = 2; %spacing
    nii = ~isnan(co(:,ii));
%    ciplot(co(nii,i)-cop(nii,i),co(nii,i)+cop(nii,i),cot(nii))
    ciplot(co(nii,ii)-cop(nii,ii) - ii*SP,co(nii,ii)+cop(nii,ii) - ii*SP,cot(nii,ii),[.7 .7 .7]);
    plot(cot(:,ii),co(:,ii) - ii*SP,'k','linewidth',2)%'color',[.7 .7 .7]);
    plot(t,ymm - ii*SP,'b.')%,'markerfacecolor','b')
   % plot(cote(:,ii),coe(:,ii),'r-','linewidth',2);

    xlim(LIMS)
%    legend('raw','interp')
    set(gcf,'color','w','position',[1076          19         560         778])
    %text(LIMS(1)-10e3,-ii*2,[core_order_new(ii);strsplit(species,',')'])
    text(-24e3,-ii*SP,[core_order_new(ii)],'fontsize',14)
    set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})

end

mvec = flipud(d18mv(:));
%set(gca,'yticklabel',{'0','0','0','0','0','0','0','0','0'});
set(gca,'yticklabel',{'';num2str(mvec,3);''});
ylabel('\delta^1^8O_c_c')
xlabel('Time (years)')

% plot([-8000 -8000],[13 14],'k')
% text(-7500,13.7,'1 permil')

gridxy(min(xlim):2500:max(xlim),min(ylim):max(ylim),'linestyle',':','color',[.5 .5 .5])

export_fig('-pdf','Figs/rawdata_inc_interp')

