function [xest] = tsvd(U,S,V,y,k)
% function [xest] = tsvd(U,S,V,y,k)
% truncated SVD solver for problems of the form y=Ex+n
%
% ---- INPUTS ----%
% 
% U, S, V are matrices from SVD of E
% y are the data
% k is the SVD truncation parameter (K' in Wunsch 2006)
%
% ---- OUTPUT ---- %
%
% xest is the estimate of the parameter x from the data (x~ in
% Wunsch 2006)
%
% D. Amrhein. Updated with comments 12 March 2015

Vk = V(:,1:k);
Uk = U(:,1:k);
Sk = S(1:k,1:k);

xest =  (Vk*(Sk\Uk'))*y(:);




