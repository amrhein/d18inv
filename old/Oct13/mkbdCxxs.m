%function [] = mkbdCxxs(k)

     k = 200;

% function [] = mkbdCxxs(k)
% a function to make and save the diagonal blocks of the whole-domain
% covariance matrix. Each of the block matrices is 2806x2806, corresponding
% to a single, spatial covariance matrix at a given time.

load 4weddell/Gbpusv23July.mat
load 4weddell/Dlsfile.mat

tau = 101;
L = 2806;
si = s(1:k,1:k)\eye(k);

Cxxs_all = nan(tau*L,L);

for ii = 1:tau
   disp(num2str(ii))
   ssi = ((ii-1)*2806+1);
  Cxxs_all(ssi,:) = v(ssi,1:k)*si*u(:,1:k)'*Pm*u(:,1:k)*si*v(ssi,1:k)';
end

fname = ['bdCxx',num2str(k)];
save(fname,'Cxxs_all')
