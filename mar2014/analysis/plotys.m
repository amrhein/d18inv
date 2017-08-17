function [] = plotys(Dwd,yestwd,Dwdt,core_order,yil,Pm)

[co] = yv2coreswd(Dwd,yil);
[coe] = yv2coreswd(yestwd,yil);
[cop] = yv2coreswd(sqrt(diag(Pm)),yil);
cote = yv2coreswd(Dwdt,yil);
cot = cote;
% cote = Dwdt;
% cot = Dwdt;
clf
set(gcf,'color','w')

% for raw data
path1 = '/home/dan/Dropbox/2012-2013/Data/benthicTMI/';
lims = [-25000,-5000];
[D,tD,Di,core_order_new] = mkDni(path1,lims);

for ii = 1:length(core_order)
    % get raw data
    ly = sum(Di==ii);
    y = D(Di==ii);
    t = tD(Di==ii);
    N = 0.2^2*eye(ly);

    ly = length(y);
    ym = get_ym(y,t,N);
    ymm = y-ym;

    
    subplot(2,4,ii)
    hold on
    nii = ~isnan(co(:,ii));
%    ciplot(co(nii,i)-cop(nii,i),co(nii,i)+cop(nii,i),cot(nii))
    ciplot(co(nii,ii)-cop(nii,ii),co(nii,ii)+cop(nii,ii),cot(nii,ii),[.7 .7 .7]);
    plot(t,ymm,'b.','markerfacecolor','b')
    plot(cot(:,ii),co(:,ii),'k','linewidth',2)%'color',[.7 .7 .7]);
    plot(cote(:,ii),coe(:,ii),'r-','linewidth',2);
%    rd = co(:,ii)-coe(:,ii);
%     plot(cote,rd,'b','linewidth',1);
%    title(char(core_order(i)),'interpreter','none','fontweight','bold')
    title([char(core_order(ii))])% , ': RMSE=' , num2str(sqrt(nanmean(rd.^2)),2)],'interpreter','none','fontweight','bold')
    xlabel('Time (years)')
    set(gca,'xtick',[-25000 -15000 -5000])
    set(gca,'xticklabel',num2str(get(gca,'xtick')',5))
    ylabel('permil \delta^1^8O_c_c')
    ylim([-1.5, 1.5])
    xlim([-25000 -5000])
    
end
set(gcf,'position',[776  -345   687   759])

addpath ../../export_fig/jlab/
%packcols(3,3)