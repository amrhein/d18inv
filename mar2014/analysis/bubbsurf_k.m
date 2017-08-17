function [s] = bubbsurf_k(vec,TOL,SCALE)
% plots just the amplitude. all colors are black.
% for times when you see a red door and you want it painted black.

addpath ../../tmi/
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
set(gcf,'color','w')%,'position',[-81         283        1727         891])
hold on

mask = s*0;
%imagescnan(lonmg(1,:),latmg(:,1),mask,'nancolor','w')
pcolor(lonmg(1,:),latmg(:,1),mask)
shading flat
colormap(0.8*ones(1,3))

niu = (~isnan(s) & s>=TOL);
a = scatter(lonmg(niu),latmg(niu),SCALE*s(niu),'k');
niun = (~isnan(s) & s<=-TOL);
aa = scatter(lonmg(niun),latmg(niun),abs(SCALE*s(niun)),'k');
axis equal, axis tight
set(a,'linewidth',0.3)
set(aa,'linewidth',0.3)

