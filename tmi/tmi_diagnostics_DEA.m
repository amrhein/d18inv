%% tmi_diagnostics.m Total Matrix Intercomparison (TMI) diagnostic routines
%  example 1: track surface waters into interior
%  example 2: Find the volume that has originated from surface
%  example 3: Find surface origin of an interior box
%  example 4: Reconstruct steady-state tracer field.
%  example 5: Construct the pathways matrix from observations.
%  example 6: Follow the adjoint transport pathways.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the master file: choose either woce2deg or woce4deg 
% (2 deg == 2 x 2 horiz resolution, 33 vertical levels)
% (4 deg == 4 x 4 horiz resolution, 33 vertical levels)
% (woce == trained on the WOCE Global Hydrographic Climatology,
%          Gouretski and Koltermann 2004)

load A_4deg_2010  %% choose a A*.mat file

% variables in A*.mat %%%
%
% A: pathways matrix
% i,j,k: vectors that contain the spatial coordinates 
%        consistent with the A matrix.
% LAT,LON,DEPTH: the latitude, longitude and depth of the 3-d grid. 
% dP = amount of remineralized phosphate in each box.

%% Most examples are computationally more efficient with the output of a LU decomposition.
% Use an LU decomposition so that the inverse does not have
% to be stored explicitly.
% ALL EXAMPLES USE THIS CALCULATION.
 [L, U, P, Q, R] = lu (A) ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1: Track the pathways of a user-defined water mass.         %
% Steps: (a) define the water mass 1). by a rectangular surface patch %
%            dyed with passive tracer concentration of 1,             % 
%            or 2. load a pre-defined surface patch in d_all.mat.     %
%        (b) propagate the dye with the matrix A, with the result     % 
%            being the fraction of water originating from the         %
%            surface region.                                          %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define the surface patch by the bounding latitude and longitude.
 lat_lo = 50; % 50 N, for example.
 lat_hi = 60;
 
 lon_lo = -50;
 lon_hi = 0;
 
 if lon_lo<=0 ;     lon_lo = lon_lo + 360; end
 if lon_hi<=0 ;     lon_hi = lon_hi + 360; end
 
 loc = find( LAT(j) < lat_hi & LAT(j)> lat_lo & ...
             LON(i) > lon_lo & LON(i)<lon_hi & k==1); 
 %potential wraparound problems in longitude are not accounted for
 %in the "loc" statement.
 
 N = length(i);
 d = zeros(N,1);
 d(loc) = 1; % a vector with ones in the N. Atl., zero elsewhere.
 
 % Also, can read in predefined surface patches in the file
 % d_all.mat, where the surface patches are defined as
 % oceanographically-relevant regions: 
 % 1) GLOBAL, 2) ANT, 3 NATL, 4 SUBANT, 5 NPAC, 
 % 6) ARC, 7 MED, 8) TROP, 9 ROSS, 10 WED, 11 LAB, 12 GIN, 13 ADEL.
 % 14) Atlantic sector SUBANT, 15) Pacific SUBANT, 16) Indian SUBANT
 % 17) Atlantic TROP, 18) Pacific TROP, 19) Indian TROP
 %
 % e.g.,
 % load d_all_4deg; d = d_all(:,3); %NATL.
 
 % to get quantity of dyed water throughout ocean:
 c = Q * (U \ (L \ (P * (R \ d))));  
 C = vector_to_field(c,i,j,k);  % g is the 3-d field of water-mass concentration.
 
 lon_section = 330;
 isec = find(LON==330);
 contourf(LAT,-DEPTH,squeeze(C(:,:,isec)),0:0.05:1) % a sample plot at 22 W.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2: Find the ocean volume that has originated from each    %
%            surface box.                                           %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of each box %
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.                                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 NZ = max(k);
 NY = max(j);
 NX = max(i);
 
 % get vector of volumes.
 dx = diff(LON(1:2));
 dy = diff(LAT(1:2));
 for ny = 1:length(LAT)
   distx(ny) = sw_dist([LAT(ny) LAT(ny)],[0 dx],'km');
 end
 disty = sw_dist([0 dy],[0 0],'km');
 distx = reshape(distx,length(distx),1); % get dimensions correct.
 area = 1e6.*disty(ones(NY,NX)).*distx(:,ones(NX,1)); 
 zface= (DEPTH(1:end-1)+DEPTH(2:end))./2;
 dz = ([zface(1) ; diff(zface); 500]);
 vol = permute(area(:,:,ones(NZ,1)),[3 1 2]) .*dz(:,ones(NY,1),ones(NX,1));
 v   = field_to_vector(vol,i,j,k);
 load inmixlyr
load tracerobs_4deg_33lev_woce.mat

%v(DEPTH(kt)>1000) = 0;
%v(inmixlyr) = 0;
 
 % effectively take inverse of transpose A matrix.
 vtot = R' \ (P' * (L' \ (U' \ (Q' * v))));  % inverse of transpose
                                             % = transpose of inverse
 
 Vtot = vector_to_field(vtot,i,j,k);
 Vtot = squeeze(Vtot(1,:,:))./area; % scale the volume by surface area.
addpath ../../export_fig/
 imagescnan(LON,LAT,log10(Vtot)) % units: effective thickness
%imagescnan(LON,LAT,(Vtot)) % units: effective thickness
set(gca,'ydir','n')
colorbar
set(gcf,'color','w')
caxis([0 7])

save Vtot_DEA_30Sept Vtot
 
