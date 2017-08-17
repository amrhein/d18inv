function [] = plotys2(Dwd,yestwd,Dwdt,core_order,yil,Pm)

addpath ../run
addpath ../setup
path1 = '../../data/benthic/';

LIMS = [-25000,-5000];
TSTEP = 200;
%load /home/dan/Dropbox/2012-2013/Data/core_order_new.mat
%load /home/dan/Dropbox/2012-2013/Data/tracerobs_4deg_33lev_woce.mat

[co] = yv2coreswd(Dwd,yil);
[coe] = yv2coreswd(yestwd,yil);
[cop] = yv2coreswd(sqrt(diag(Pm)),yil);
cote = yv2coreswd(Dwdt,yil);
cot = cote;

clf
set(gcf,'color','w')

% for raw data
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
    d18mv(ii) = ym;
    ymm = y-ym;

    hold on

SP = 2; %spacing
    nii = ~isnan(co(:,ii));
    ciplot(co(nii,ii)-cop(nii,ii) - ii*SP,co(nii,ii)+cop(nii,ii) - ii*SP,cot(nii,ii),[.7 .7 .7]);
    plot(cot(:,ii),co(:,ii) - ii*SP,'k','linewidth',2)%'color',[.7 .7 .7]);
    plot(t,ymm - ii*SP,'b.')%,'markerfacecolor','b')
    plot(cote(:,ii),coe(:,ii) - ii*SP,'r-','linewidth',2);

    xlim(LIMS)

    set(gcf,'color','w','position',[748    62   337   752])
    %    text(-24e3,-ii*SP,[core_order_new(ii)],'fontsize',14)
    text(-24e3,-ii*SP,[core_order(ii)],'fontsize',14,'interpreter','none')
    set(gca,'xticklabel',{'-25000','-20000','-15000','-10000','-5000'})

end

mvec = flipud(d18mv(:));
set(gca,'yticklabel',{'0','0','0','0','0','0','0','0','0'});
%set(gca,'yticklabel',{'';num2str(mvec,3);''});
ylabel('\delta^1^8O_c anomalies (permil PDB)')
xlabel('Time (years)')


gridxy(min(xlim):2500:max(xlim),min(ylim):max(ylim),'linestyle',':','color',[.5 .5 .5])



