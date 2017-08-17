%function [] = mkbdCxxs(k)

     k = 200;

% function [] = mkbdCxxs(k)
% a function to make and save the diagonal blocks of the whole-domain
% covariance matrix. Each of the block matrices is 2806x2806, corresponding
% to a single, spatial covariance matrix at a given time.

load 4weddell/Gbpusv23July.mat
load 4weddell/Dlsfile.mat
load c_all_4deg.mat  
load Rig.mat  
load tracerobs_4deg_33lev_woce.mat  
load Vtot_DEA_30Sept.mat

Vtotv = Vtot(~isnan(Vtot));
L = length(Vtotv);

si = s\eye(size(s));

c_all_surf = c_all(kt==1,:)';

nr = min(size(c_all_surf));

tau = 101;
L = 2806;
kmax = [200];

Crr = nan(tau,kmax,nr);

tr = trace(c_all_surf*sparse(diag(Vtotv)));

for ii = 1:tau
   ti = ((ii-1)*L+1);
   tti = ti:(ti+L-1);

   for k = 1:kmax

      % Compute the covariance contribution from a single singular vector index over
      % the time specified by ssi
      % Cxx = v(ssi,k)*inv(s(k,k))*u(ssi,k)'*Pm*u(ssi,k)'*inv(s(k,k))*v(ssi,k)'
      
      % Compute the contribution for each region
      % HCxxH = sparse(diag(Vtotv))*...
      %   v(ssi,k)*inv(s(k,k))*u(ssi,k)'*Pm*u(ssi,k)'*inv(s(k,k))*v(ssi,k)'...
      %   *sparse(diag(Vtotv));

      Crr(ii,k,:) = tr^-2*c_all_surf*...
         sparse(diag(Vtotv))*...
         v(tti,k)*si(:,k)*u(:,k)'*Pm*u(:,k)'*si(:,k)*v(tti,k)'...
         *sparse(diag(Vtotv))...
         *c_all_surf';

   end
end
