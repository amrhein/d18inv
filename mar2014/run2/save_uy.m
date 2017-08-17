load Gbpusv_sqrt u s
load ../setup/Dlsfile

kk = 806;
kkv = 1:kk;

sd=diag(s(1:kk,1:kk));

% projection of data onto singular vectors
uy = (u(:,1:kk)'*Dwd);
% soln variance
sv = cumsum((inv(s(1:kk,1:kk))*u(:,1:kk)'*Dwd).^2);
% amt by which sum of squared residuals is reduced
ssrr = sum(Dwd.^2)-cumsum((u(:,1:kk)'*Dwd).^2);

% see 2.312 Wunsch 2006
%ns = bsxfun(@times,(u'*inv(chol(Pm))*Dwd(:))',u);
ns = bsxfun(@times,(u'*Dwd(:))',u);
nc = fliplr(cumsum(fliplr(ns),2));
J = sum(nc.^2);
Js = sum(((chol(Pm))\nc).^2);

save uy.mat kk kkv sd uy sv ssrr
