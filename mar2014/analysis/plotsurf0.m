function [map] = plotsurf0(vec)
% makes the 0 values black

addpath ../../tmi/
addpath ../../../export_fig/
load tracerobs_4deg_33lev_woce it jt kt

if length(vec) == length(kt)
    vec1 = vec;
elseif length(vec) == sum(kt==1)
    vec1 = nan(size(kt));
    vec1(kt==1) = vec;
else
    error('bad vector length')
end

vf = vector_to_field(vec1,it,jt,kt);
map = squeeze(vf(1,:,:));
LAT = -88:4:88;
LON = 2:4:358;
clf
imagescnan(LON,LAT,map,'nancolor','w');
%pcolor(LON,LAT,map)
set(gca,'ydir','n')
cax = caxis;
%cm0 = redblue(32);
CR= 20; % fix to make sure that the lowest values aren't zeroed out
cm0 = b2r([cax(1)-diff(cax)/CR,cax(2)],64);
map0 = map;
%map0(map==0)=-9999;
map0(abs(map)<eps)=-9999;
%imagescnan(LON,LAT,map0,'nancolor','w')
imagescnan(LON,LAT,map0,'nancolor',[0.5 0.5 0.5])
colormap([[0.7 0.7 0.7];cm0])
caxis([cax(1)-diff(cax)/CR,cax(2)])
colorbar
%colorbar('westoutside')
%shading flat
set(gca,'ydir','n')

