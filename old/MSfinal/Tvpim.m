% a script to filter the planktonic records using Tv. Tv will not be stored in memory so as to save time.

clear

load ../MSfinal/4weddell/Gbpusv.mat
load pim

kk = 413;
NS = 2806; % number of surface points

[r, c] = size(pim);

for ii = 1:c
   ind = tmi_ind(ii);
   pi = pim(:,ii);
   pnn = ~isnan(pi);
   p = pi(pnn)
   vs = v(ind:NS:end,1:kk); % subselected in space
   vst = vs(pnn,:); % subselected in time
   vp = vst'*p;
   tvp(pnn,ii) = vst*vp;
end

save tvp tvp
