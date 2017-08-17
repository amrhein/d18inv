
load ../MSfinal/4weddell/Gbpusv.mat u s
load ../MSfinal/4weddell/Dlsfile.mat

% clf
% set(gcf,'color','w','position',[1318         170        1027         418])
% kk = length(Dwd);
% kkv = 1:kk;
% 
% plot(kkv,diag(s(1:kk,1:kk)))
% hl = legend('\lambda_i');
% set(gca,'fontsize',16)
% ylabel('Normalized units','fontsize',16)
% xlabel('SVD index','fontsize',16)
% %print -dpdf -r300 Figs/uy1
% 
% clf
% [ax,h1,h2]=plotyy(kkv,abs(u(:,1:kk)'*Dwd),kkv,abs(inv(s(1:kk,1:kk))*u(:,1:kk)'*Dwd))
% hold on
% plot(xlim,[1 1],':','color',[0.5 0.5 0.5])
% hl = legend('u_i^Ty','u_i^Ty/\lambda_i');
% set(hl,'location','northwest')%,'box','off')
%      ylim(ax(1),[0,5])
%      ylim(ax(2),[0,5])
%      xlabel(ax(1),'SVD index','fontsize',16)
%      xlabel(ax(2),'SVD index','fontsize',16)
%      ylabel(ax(1),'Concentration','fontsize',16)
%      ylabel(ax(2),'Normalized units','fontsize',16)
% set(ax(1),'fontsize',16)
% set(ax(2),'fontsize',16)

%print -dpdf -r300 Figs/uy2

clf
set(gcf,'position',[1318           5         446         792],'color','w')
%% Variance of the solution
subplot(4,1,3)
plot(kkv,cumsum((inv(s(1:kk,1:kk))*u(:,1:kk)'*Dwd).^2))
ylim([0 5e3])
xlabel('SVD index')
ylabel('Concentration squared')
xlim([0 785])
ylabel(['Solution variance (' char(8240) ')'])
%% amt by which sum of squared residuals is reduced
subplot(4,1,4)
plot(kkv,sum(Dwd.^2)-cumsum((u(:,1:kk)'*Dwd).^2))
ylim([0, sum(Dwd.^2)])
xlim([0 785])
xlabel('SVD index')
ylabel('Residual variance')
%% projection onto singular vectors
subplot(4,1,1)
plot(kkv,(u(:,1:kk)'*Dwd))
xlim([0 785])
xlabel('SVD index')
ylabel('u_i^Ty')
%%
subplot(4,1,2)
plot(kkv,diag(inv(s(1:kk,1:kk))).^2)
ylim([0 1000])
xlim([0 785])
xlabel('SVD index')
ylabel('1/\lambda_i')
%%


