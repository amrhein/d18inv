function [map] = plotsurf(vec)

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
imagescnan(LON,LAT,map,'nancolor','w')
%pcolor(LON,LAT,map)
set(gca,'ydir','n')
%colorbar('westoutside')
%shading flat
    
