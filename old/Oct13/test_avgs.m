%% setup
%load 4weddell/Gbpusv23July.mat
%load Vtot_DEA_30Sept          
%load c_all_4deg 

%%
tau  = 101;

Vtotv = Vtot(~isnan(Vtot));
L = length(Vtotv);
ii = 1;
surfi = plotsurf(c_all(:,ii));
surfis = surfi(~isnan(surfi(:))).*Vtotv;

%Ri = sptoeplitz([surfis(1);zeros(tau-1,1)],[surfis;zeros((tau-1)*L,1)]);

Ri = sparse(tau,tau*L);
for jj = 1:tau
   si = L*(jj-1);
Ri(jj,(si+1):(si+L)) = surfis;
end

save Rig Ri 
