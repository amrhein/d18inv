
%load /home/dan/Dropbox/2012-2013/MSfinal/4weddell/Gbpusv.mat 
%load /home/dan/Dropbox/2012-2013/MSfinal/4weddell/Dlsfile.mat
clear
load fweddell/icout3.mat
kk = 8;
kkv = 1:kk;
uy = (uic'*yold);
sv = cumsum((inv(sic(1:kk,1:kk))*uic'*yold).^2);
ssrr = sum(yold.^2)-cumsum((uic'*yold).^2);
sd = diag(sic);

figure(1)
clf
set(gcf,'position',[1318           5         446         792],'color','w')

%% yet another way

figure(2)
clf
[ax,h1,h2] = plotyy(kkv,sv,kkv,ssrr);
xlabel(ax(1),'$K^\prime$','interpreter','latex')
set(gca,'position',[ 0.1300    0.1700    0.7750    0.8050])
grid on
%set(ax(2),'xlim',[0 785],'ylim',[0,400],'ytick',linspace(0,400,5))
%set(ax(1),'xlim',[0 785],'ylim',[0,5e3],'ytick',linspace(0,5e3,5))
% set(h1,'linewidth',2)
% set(h2,'linewidth',2)
set(gcf,'color','w','position',[1177         335         560         229])
hold on
plot([3 3],ylim(gca),'--k')
ylabel(ax(1),'Solution RMS (permil)')
ylabel(ax(2),'Residual RMS (permil)')
export_fig('-pdf','Figs/uy1_ic')

figure(3)
clf
[ax,h1,h2] = plotyy(kkv,uy,kkv,1./sd);
xlabel(ax(1),'SVD index')
set(gca,'position',[ 0.1300    0.1700    0.7750    0.8050])
%set(ax(1),'xlim',[0 785],'ylim',[-7 7],'ytick',-6:3:6)
%set(ax(2),'xlim',[0 785],'ylim',[0,1e2],'ytick',0:20:1e2)
set(gcf,'color','w','position',[1177         335         560         229])
hold on
plot([3 3],ylim(gca),'--k')
grid
ylabel(ax(1),'Projection of y onto u_i (permil)')
ylabel(ax(2),'1/\lambda_i')
export_fig('-pdf','Figs/uy2_ic')


% % projection onto singular vectors
% subplot(4,1,1)
% plot(kkv,(uic(:,1:kk)'*yold))
% %xlim([0 785])
% xlabel('SVD index')
% ylabel('u_i^Ty')
% %
% subplot(4,1,2)
% plot(kkv,diag(inv(sic(1:kk,1:kk))).^2)
% %ylim([0 1000])
% %xlim([0 785])
% xlabel('SVD index')
% ylabel('1/\lambda_i')
% 
% % Variance of the solution
% kk = 8;
% kkv = 1:kk;
% subplot(4,1,3)
% plot(kkv,cumsum((inv(sic(1:kk,1:kk))*uic(:,1:kk)'*yold).^2))
% xlabel('SVD index')
% ylabel('Concentration squared')
% 
% % amt by which sum of squared residuals is reduced
% subplot(4,1,4)
% plot(kkv,sum(yold.^2)-cumsum((uic(:,1:kk)'*yold).^2))
% %ylim([0, sum(yic.^2)])
% %xlim([0 785])
% xlabel('SVD index')
% ylabel('Residual variance')
% 
% 
% 
% export_fig('-pdf',['Figs/xic_svd_inds',num2str(kic)])


