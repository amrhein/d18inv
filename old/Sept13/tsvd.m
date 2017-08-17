function [xest] = tsvd(U,S,V,y,k)
% function [xest] = tsvd(U,S,V,y,sig,k)
% truncated SVD solver

Vk = V(:,1:k);
Uk = U(:,1:k);
Sk = S(1:k,1:k);

xest =  (Vk*(Sk\Uk'))*y(:);

% % compute the diagonals of the solution covariance matrix
% Vks = Vk*inv(Sk.^2);
% dCxx = nan(length(Vk),1);
% matlabpool open 4
% parfor ii = 1:length(Vk)
% for ii = 1:length(Vk)
%     dCxx(ii) = (Vks(ii,:)*Vk(ii,:)');
% end
% dCxx = dCxx*sig^2;
%dCxx = sig^2*Vk*inv(Sk.^2)*Vk';
%tic,E1 = sum((T*H).*T,2);toc

%Tu = Uk*Uk';
%Tv = Vk*Vk';



