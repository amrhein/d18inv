%function [] = mkbdCxxs(k)
% function [] = mkbdCxxs(k)
% a function to make and save the diagonal blocks of the whole-domain
% covariance matrix. Each of the block matrices is 2806x2806, corresponding
% to a single, spatial covariance matrix at a given time.

clear

k = 400;

load 4weddell/Gbpusv23July.mat
load 4weddell/Dlsfile.mat
load c_all_4deg.mat  
load tracerobs_4deg_33lev_woce.mat  
load Vtot_DEA_30Sept.mat

Vtotv = Vtot(~isnan(Vtot));
L = length(Vtotv);

si = s\eye(size(s));

c_all_s = c_all(kt==1,:)';
c_all_surf = c_all_s(1:2,:);

nr = min(size(c_all_surf));

tau = 101;
L = 2806;

Crr = nan(tau,nr);

tr = c_all_surf*sparse(diag(Vtotv));

for ii = 1:tau
disp(['t = ' num2str(ii)])
   ti = ((ii-1)*L+1);
   tti = ti:(ti+L-1);

Crr(ii,:) = diag(...
         c_all_surf*...
         sparse(diag(Vtotv))*...
         v(tti,1:k)*si(1:k,1:k)*(u(:,1:k)'*Pm*u(:,1:k))*si(1:k,1:k)*v(tti,1:k)'...
         *sparse(diag(Vtotv))...
         *c_all_surf');

end

Crr = Crr./repmat(sum(tr').^2,tau,1);
