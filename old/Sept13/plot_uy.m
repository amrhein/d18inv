load ../23july/4weddell/Gbpusv23July.mat u s
load ../23july/4weddell/Dlsfile.mat

clf
set(gcf,'color','w')
%hold all

kk = 500;
kkv = 1:kk;

[ax,h1,h2] = plotyy(kkv,diag(s(1:kk,1:kk)),kkv,inv(s(1:kk,1:kk))*u(:,1:kk)'*Dwd);
hold(ax(1),'on');hold(ax(2),'on');
h3 = plot(ax(2),kkv,u(:,1:kk)'*Dwd,'r');
%plot(diag(s(1:kk,1:kk)))
hl = legend([h1,h2,h3],'u_i^Ty/\lambda_i','\lambda_i','u_i^Ty');
set(hl,'box','off')



