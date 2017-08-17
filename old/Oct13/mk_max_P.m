% mk_max_P.m
% a script to determine, using svd, the maximum projection vector v
% onto the regional-averaged concentration output for a given region.  

clear

k = 400;

load 4weddell/Gbpusv23July.mat
load 4weddell/Dlsfile.mat
load c_all_4deg.mat  
load Rig
load tracerobs_4deg_33lev_woce.mat  
load Vtot_DEA_30Sept.mat

Vtotv = Vtot(~isnan(Vtot));
L = length(Vtotv);


c_all_s = c_all(kt==1,:)';
c_all_surf = c_all_s(1,:);

nr = min(size(c_all_surf));

tau = 101;
L = 2806;

Crr = nan(tau,nr);

tr = c_all_surf*sparse(diag(Vtotv));


%App = c_all_surf*sparse(diag(Vtotv))- ...
%         (c_all_surf*...
%         sparse(diag(Vtotv))*...
%         v(:,1:k))*v(:,1:k)';
         
%App = Ri-(Ri*v(:,1:k))*v(:,1:k)';

% see what happens if I eliminate the nullspace at the end of the recon
%Ri(:,(end-5*L+1):end)=[];
%v((end-5*L+1):end,:)=[];  


[ur sr vr] = svd((Ri-(Ri*v(:,1:k))*v(:,1:k)')/sum(tr),'econ');

%[ujr sjr vjr] = svd(full(Ri)/sum(tr),'econ');
%vrr = [reshape(vr,2806,[])]; 
%[u1,s1,v1] = svd(vrr);

%Crr = Crr./repmat(sum(tr').^2,tau,1);
