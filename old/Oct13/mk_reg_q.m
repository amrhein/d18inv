 function [r] = mk_reg_q(q)
% function [regq] = mk_reg_q(q)
% a function to make and save the diagonal blocks of the whole-domain
% covariance matrix. Each of the block matrices is 2806x2806, corresponding
% to a single, spatial covariance matrix at a given time.

load c_all_4deg.mat  
load tracerobs_4deg_33lev_woce.mat  
load Vtot_DEA_30Sept.mat

Vtotv = Vtot(~isnan(Vtot));
L = length(Vtotv);

c_all_s = c_all(kt==1,:)';

c_all_surf = c_all_s;

nr = min(size(c_all_surf));

tau = 101;
L = 2806;
kmax = [200];

tr = c_all_surf*sparse(diag(Vtotv));

for ii = 1:tau
   disp(['t = ' num2str(ii)])

   ti = ((ii-1)*L+1);
   tti = ti:(ti+L-1);

   r(ii,:) = tr*q(tti);
   
end

r = r./repmat(sum(tr'),tau,1);
