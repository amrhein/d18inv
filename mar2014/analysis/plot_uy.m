load uy.mat
figure(1)
clf
set(gcf,'position',[1318           5         446         792],'color','w')
%% Variance of the solution
subplot(4,1,3)
plot(kkv,sv)
ylim([0 5e3])
xlabel('SVD index')
xlim([0 785])
ylabel(['Solution RMS (permil)'])
%% amt by which sum of squared residuals is reduced
subplot(4,1,4)
plot(kkv,ssrr)
%ylim([0, sum(Dwd.^2)])
xlim([0 785])
xlabel('SVD index')
ylabel('Residual variance')
%% projection onto singular vectors
subplot(4,1,1)
plot(kkv,uy)
xlim([0 785])
xlabel('SVD index')
ylabel('u_i^Ty')
axis tight
%%
subplot(4,1,2)
plot(kkv,1./sd)
ylim([0 100])
xlim([0 785])
xlabel('SVD index')
ylabel('1/\lambda_i')

%%
addpath ../../../export_fig/
export_fig('-pdf','Figs/uy')


%% another way
figure(2)
clf
[ax,hlines] = plotyyy(kkv,ssrr,kkv,1./sd,kkv,sv);
set(ax(1),'xlim',[0 785],'ylim',[0,400],'ytick',linspace(0,400,9))
set(ax(2),'xlim',[0 785],'ylim',[0,1e2],'ytick',0:10:1e2)
set(ax(3),'xlim',[0 785],'ylim',[0,5e3])

%% yet another way

figure(2)
clf
[ax,h1,h2] = plotyy(kkv,sv,kkv,ssrr);
xlabel(ax(1),'$K^\prime$','interpreter','latex')
set(gca,'position',[ 0.1300    0.1700    0.7750    0.8050])
grid on
set(ax(2),'xlim',[0 785],'ylim',[0,400],'ytick',linspace(0,400,5))
set(ax(1),'xlim',[0 785],'ylim',[0,5e3],'ytick',linspace(0,5e3,5))
% set(h1,'linewidth',2)
% set(h2,'linewidth',2)
set(gcf,'color','w','position',[1177         335         560         229])
hold on
plot([413 413],ylim(gca),'--k')
ylabel(ax(1),'Solution RMS (permil)')
ylabel(ax(2),'Residual RMS (permil)')
export_fig('-pdf','Figs/uy1')

figure(3)
clf
[ax,h1,h2] = plotyy(kkv,uy,kkv,1./sd);
xlabel(ax(1),'SVD index')
set(gca,'position',[ 0.1300    0.1700    0.7750    0.8050])
set(ax(1),'xlim',[0 785],'ylim',[-7 7],'ytick',-6:3:6)
set(ax(2),'xlim',[0 785],'ylim',[0,1e2],'ytick',0:20:1e2)
set(gcf,'color','w','position',[1177         335         560         229])
hold on
plot([413 413],ylim(gca),'--k')
grid
ylabel(ax(1),'Projection of y onto u_i (permil)')
ylabel(ax(2),'1/\lambda_i')
export_fig('-pdf','Figs/uy2')
