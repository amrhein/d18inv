% mkG.m
% This function constructs the observational design matrix (here called G)
% from the Green's function matrix g for a specified time stlims
% and time step TSTEP.

function [G] = mkG_2(stlims,dtlims,TSTEP,g,m,l)

% g has size [tau,n] where tau is the number of times and n = m*l is the
% number of observations times the number of controls. the columns of g are
% ordered first by source region and then by observational locus.

% construct the first block column of G. The number of block rows is equal
% to the number of time steps. I need to append a zero block matrix at the
% top corresponding to the fact that no interior locations see the boundary
% conditions instantaneously.

nt = abs(diff(stlims)/TSTEP)+1;
ndt = abs(diff(dtlims)/TSTEP)+1;
G0 = zeros(m,l);
G1 = G0;

[ngt,blah] = size(g);

% first create the part of the block column with block length equal to the
% number of time steps in the Green's functions
for ii = 1:ngt
    G1 = [G1;reshape(g(ii,:),l,m)'];
end

% now adjust the block column length according to stlims by adding zeros
if     nt < ngt
    G1(((nt-1)+1)*m:end,:) = [];
elseif nt > ngt
    G1 = [G1;repmat(G0,(nt-ngt-1),1)];
end
   
% Now row-concatenate additional block-column vectors

G = G1;
for ii = 1:dt-1
    G = [G,[repmat(G0,ii,1);G1(1:(end-ii*m),:)]];
end




