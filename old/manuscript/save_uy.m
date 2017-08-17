load ../MSfinal/4weddell/Gbpusv.mat u s
load ../MSfinal/4weddell/Dlsfile.mat

kk = 785;
kkv = 1:kk;

sd=diag(s(1:kk,1:kk));

% projection of data onto singular vectors
uy = (u(:,1:kk)'*Dwd);
% soln variance
sv = cumsum((inv(s(1:kk,1:kk))*u(:,1:kk)'*Dwd).^2);
% amt by which sum of squared residuals is reduced
ssrr = sum(Dwd.^2)-cumsum((u(:,1:kk)'*Dwd).^2);

save uy.mat kk kkv sd uy sv ssrr
