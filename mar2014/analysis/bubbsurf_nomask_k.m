function [s] = bubbsurf_nomask_k(vec,TOL,SCALE)

addpath ../../tmi
load ../../tmi/tracerobs_4deg_33lev_woce it jt kt

if length(vec) == length(kt)
    vec1 = vec;
elseif length(vec) == sum(kt==1)
    vec1 = nan(size(kt));
    vec1(kt==1) = vec;
else
    error('bad vector length')
end

vf = vector_to_field(vec1,it,jt,kt);
s = squeeze(vf(1,:,:));

[lonmg,latmg] = meshgrid(2:4:358,-90:4:88);

%close all
%set(gcf,'color','w')%,'position',[-81         283        1727         891])
%hold on

mask = s*0;
%imagescnan(lonmg(1,:),latmg(:,1),mask) %mask
%colormap bone

niu = (~isnan(s) & s>=TOL);
scatter(lonmg(niu),latmg(niu),SCALE*s(niu),'k')
niun = (~isnan(s) & s<=-TOL);
scatter(lonmg(niun),latmg(niun),abs(SCALE*s(niun)),'k')
axis equal, axis tight

