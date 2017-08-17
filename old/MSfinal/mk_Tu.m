function tvi = mk_Tu(V,k,ii)
% function tui = mk_Tu(U,k,i)
% a function to return one row (or column) of the symmetric resolution
% matrix Tu = Uk*Uk'. U is the matrix of n singular vectors of length m
%(arranged columnwise), k<min(m,n) is the nullspace parameter, and i
%specifies the row/column of interest.

%keyboard
tvi = V(:,1:k)*V(ii,1:k)';

end